import pathlib
import sys
module_dir = pathlib.Path(__file__).parent.resolve()
root_dir = module_dir.parent
asvlite_wrapper_dir = root_dir.joinpath("dependency", "ASVLite", "wrapper", "cython")
sys.path.insert(0, str(asvlite_wrapper_dir))
import math
import multiprocessing as mp
import pyproj
import epsg
import random
import numpy as np
from tqdm import tqdm # to show a progress bar
from sea_surface import py_Sea_surface
from asv import py_Asv, py_Asv_specification
from geometry import py_Coordinates_3D

class Rudder_PID_controller:
    def __init__(self, asv_spec, K=None): 
        self.max_rudder_angle = math.pi/6 # 30 deg
        self.asv_spec = asv_spec
        self.error = 0.0 # error in each time step
        self.previous_error  = 0.0 # error in the previous time step
        self.cumulative_error = 0.0 # integral of errors
        self.delta_error = 0.0 # differential of error
        self.out_dir = root_dir.joinpath("results", "rudder_controller", "tuning")
        self.out_dir.mkdir(parents=True, exist_ok=True) # Create the out_dir if it does not exist.
        # P,I,D gain terms
        if K == None:
            initial_PIDs = []
            for P in [random.randint(1,3), random.randint(1,3)]:
                for I in [random.randint(1,3), random.randint(1,3)]:
                    for D in [random.randint(1,3), random.randint(1,3)]:
                        initial_PIDs.append([P,I,D])
            initial_PIDs = tqdm(initial_PIDs, leave=False)   
            initial_PIDs.set_description("Tuning iterations")  
            for initial_PID in initial_PIDs:
                P,I,D = initial_PID
                self.K = np.array([P,I,D])
                self._tune_controller()
        else:
            self.K = np.array(K)

    def __relative_angle(self, asv, waypoint):
        theta = None
        p1 = asv.py_get_position_origin()
        p2 = asv.py_get_position_cog()
        p3 = waypoint
        # theta = math.atan2((p3.y-p1.y), (p3.x-p1.x)) - math.atan2((p2.y-p1.y), (p2.x-p1.x))
        long1, lat1 = epsg.PCS_to_GCS(p1.x, p1.y)
        long2, lat2 = epsg.PCS_to_GCS(p2.x, p2.y)
        long3, lat3 = epsg.PCS_to_GCS(p3.x, p3.y)
        geodesic = pyproj.Geod(ellps='WGS84')
        fwd_azimuth_1, back_azimuth_1, distance_1 = geodesic.inv(long1, lat1, long2, lat2)
        # fwd_azimuth_1 = fwd_azimuth_1 if fwd_azimuth_1 >= 0.0 else (360 + fwd_azimuth_1)
        fwd_azimuth_2, back_azimuth_2, distance_2 = geodesic.inv(long1, lat1, long3, lat3)
        # fwd_azimuth_2 = fwd_azimuth_2 if fwd_azimuth_2 >= 0.0 else (360 + fwd_azimuth_2)
        theta = (fwd_azimuth_2 - fwd_azimuth_1) * math.pi/180
        return theta
       
    def get_rudder_angle(self, asv, waypoint):
        # Compute the relative angle between the vehicle heading and the waypoint.
        theta = self.__relative_angle(asv, waypoint)
        # Set error as the difference of the current heading and the desired heading.
        self.previous_error = self.error
        self.error = theta
        gamma = 0.7 # Rate at which the past errors reduces. 
        self.cumulative_error = self.error + gamma*self.cumulative_error
        self.delta_error = self.error - self.previous_error
        # Compute the rudder angle
        E = np.array([self.error, self.cumulative_error, self.delta_error]) # P, I, D errors.
        phi = np.dot(np.transpose(self.K), E) # radians because error is in radians.
        # Limit the rudder angle within the range (-PI/6, PI/6)
        if phi > self.max_rudder_angle:
            phi = self.max_rudder_angle
        elif phi < -self.max_rudder_angle:
            phi = -self.max_rudder_angle
        return phi # radians
   
    def _simulate_asv_in_sea_state(self, sea_state_and_PID):
        significant_wave_ht, asv_heading, P, I, D = sea_state_and_PID
        long, lat = (-150.0, 20)
        x, y = epsg.GCS_to_PCS(long, lat)
        start_point = py_Coordinates_3D(x, y, 0)
        long, lat = (-150.0, 20.01)
        x, y = epsg.GCS_to_PCS(long, lat)
        waypoint = py_Coordinates_3D(x, y, 0)
        # Init waves
        rand_seed = 1
        count_wave_spectral_directions = 5
        count_wave_spectral_frequencies = 7
        wave = py_Sea_surface(significant_wave_ht, 0.0, rand_seed, count_wave_spectral_directions, count_wave_spectral_frequencies)
        # Init ASV
        time_step_size = 40 # milli sec
        attitude = py_Coordinates_3D(0.0, 0.0, asv_heading)
        asv = py_Asv(self.asv_spec, wave, start_point, attitude)
        thrust_tuning_factor = 0.03 # The thrust tuning factor is an assumed and, for controller tuning, thrust tuning factor 
                                    # is assumed as a constant for all sea states. 
        asv.py_wg_set_thrust_tuning_factor(thrust_tuning_factor) 
        # Init controller
        controller = Rudder_PID_controller(self.asv_spec, [P,I,D])
        # Simulate
        time = 0.0
        max_time = 60 # sec
        heading_error = 0
        while time < max_time:
            time += time_step_size/1000.0
            # Compute the dynamics
            rudder_angle = controller.get_rudder_angle(asv, waypoint)
            asv.py_wg_compute_dynamics(rudder_angle, time_step_size)
            # Compute the error in heading
            error = self.__relative_angle(asv, waypoint)
            heading_error += error*error
        # Root mean square error
        rms_error = math.sqrt(heading_error/(max_time * (1000/time_step_size)))
        return rms_error
    
    def _tune_controller(self):  
        f = open("{}/{}_{}_{}.txt".format(self.out_dir, self.K[0], self.K[1], self.K[2]), "w")   
        f.write("P I D error_avg error_std\n")
        pool = mp.Pool(mp.cpu_count()) # Create a pool of processes to run in parallel. 
        delta = 0.25
        P_current, I_current, D_current = list(self.K)     
        # Policy Gradient Reinforcement Learning
        iterations = tqdm(range(10), leave=False) # This is going to take some time, therefore show a progress bar.
        iterations.set_description("Policy iterations")
        for n in iterations: 
            costs = []
            PIDs = []
            for P in [P_current-delta, P_current, P_current+delta]:
                for I in [I_current-delta, I_current, I_current+delta]:
                    for D in [D_current-delta, D_current, D_current+delta]:
                        PIDs.append([P,I,D])
            PIDs = tqdm(PIDs, leave=False)   
            PIDs.set_description("Controller variants")         
            for PID in PIDs:
                P,I,D = PID
                sea_states_and_PID = []
                for significant_wave_ht in np.arange(1.0, 10.0, 2.0):
                    for asv_heading in np.arange(0.0, 360.0, 45.0):
                        sea_states_and_PID.append([significant_wave_ht, asv_heading * math.pi/180, P, I, D])
                results = [] 
                for result in pool.imap_unordered(self._simulate_asv_in_sea_state, sea_states_and_PID): # Run multiple simulations in parallel
                    results.append(result) # append the return for each call to self._simulate_asv_in_sea_state to the list. 
                # Compute cost for each combination of PID:
                costs.append([P, I, D, np.average(np.array(results))])
            # Compute the next set of PID terms
            costs = np.array(costs)
            f.write("{} {} {} {} {}\n".format(  P_current, 
                                                I_current, 
                                                D_current, 
                                                np.average(costs, axis=0)[-1], 
                                                np.std(costs, axis=0)[-1]))
            def compute_average_cost(index, K):
                mask = []
                for item in costs:
                    mask_value = True if item[index] == K else False
                    mask.append(mask_value)
                average_cost = costs[mask].mean(axis=0)[-1]
                return average_cost
            # Compute the average performance for all cases with P_current-delta
            avg_costs_P_minus = compute_average_cost(0, P_current-delta)
            # Compute the average performance for all cases with P_current
            avg_costs_P = compute_average_cost(0, P_current)
            # Compute the average performance for all cases with P_current+delta
            avg_costs_P_plus = compute_average_cost(0, P_current+delta)
            # Compute the average performance for all cases with I_current-delta
            avg_costs_I_minus = compute_average_cost(1, I_current-delta)
            # Compute the average performance for all cases with I_current
            avg_costs_I = compute_average_cost(1, I_current)
            # Compute the average performance for all cases with I_current+delta
            avg_costs_I_plus = compute_average_cost(1, I_current+delta)
            # Compute the average performance for all cases with D_current-delta
            avg_costs_D_minus = compute_average_cost(2, D_current-delta)
            # Compute the average performance for all cases with D_current
            avg_costs_D = compute_average_cost(2, D_current)
            # Compute the average performance for all cases with D_current+delta
            avg_costs_D_plus = compute_average_cost(2, D_current+delta)
            # Compute the Adjustment vector.
            A = [0,0,0] 
            # Adjustment for P
            min_costs_P = min(avg_costs_P_minus, avg_costs_P, avg_costs_P_plus)
            if min_costs_P == avg_costs_P_minus:
                A[0] = -1
            elif min_costs_P == avg_costs_P:
                A[0] = 0
            else:
                A[0] = 1
            # Adjustment for I
            min_costs_I = min(avg_costs_I_minus, avg_costs_I, avg_costs_I_plus)
            if min_costs_I == avg_costs_I_minus:
                A[1] = -1
            elif min_costs_I == avg_costs_I:
                A[1] = 0
            else:
                A[1] = 1
            # Adjustment for D
            min_costs_D = min(avg_costs_D_minus, avg_costs_D, avg_costs_D_plus)
            if min_costs_D == avg_costs_D_minus:
                A[2] = -1
            elif min_costs_D == avg_costs_D:
                A[2] = 0
            else:
                A[2] = 1
            # Compute the new gain terms
            A = np.array(A)
            K_current = np.array([P_current, I_current, D_current])
            K_current = K_current + A*delta
            P_current, I_current, D_current = list(K_current)
            self.K = K_current
        f.close()

if __name__ == '__main__': 
    import cProfile
    # Wave glider specs
    asv_spec = py_Asv_specification()
    asv_spec.L_wl       = 2.1 # m
    asv_spec.B_wl       = 0.6 # m
    asv_spec.D          = 0.25 # m 
    asv_spec.T          = 0.09 # m
    asv_spec.max_speed  = 4.0 # m/s
    asv_spec.disp       = 0.09 # m3
    asv_spec.r_roll     = 0.2 # m
    asv_spec.r_pitch    = 0.6 # m
    asv_spec.r_yaw      = 0.6 # m
    asv_spec.cog        = py_Coordinates_3D(1.05, 0.0, -3.0) # m
    rudder_controller   = Rudder_PID_controller(asv_spec) # Will also tune the controller. 