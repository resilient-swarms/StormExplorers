import os
import re
import requests
import pathlib
import pandas as pd
from datetime import datetime
from bs4 import BeautifulSoup
from datetime import datetime

class Storms_archieve:
    '''
    This class provides interface to access the NOAA's web archive to create a local 
    database of tropical cyclones in the North Atlantic.
    '''
    def __init__(self, file_path, ignore_local_file=False):
        '''
        Initialise a local archive of storms if it does not exist.  
        '''
        self.file_path = file_path
        local_db_exists = None
        if os.path.exists(file_path):
            local_db_exists = True 
        else:
            local_db_exists = False
        # Gather data from online archive only if necessary.
        if local_db_exists == False or ignore_local_file == True:
            f = open(file_path, "w")
            f.write("year storm_name time_stamps path(long,lat) min_atm_pressure(MB) max_sustained_windspeed(knots) predicted_paths(long,lat)\n".replace(" ", "\t"))    
            # Access NOAA web archive and create a local database of North Atlantic storms
            for year in range(1999, datetime.now().year+1, 1):
                url_for_each_year = "https://www.nhc.noaa.gov/archive/{}/".format(year)
                archive_html = self.__get_html(url_for_each_year.format(year)) if year != 2002 else self.__get_html("https://www.nhc.noaa.gov/archive/2002/index.shtml")
                archive_soup = BeautifulSoup(archive_html, "html.parser")
                # Html page contains a table with a list of tropical storms for the year. 
                # In general this will be the only table in the page, but for the years 1999, 2000, 2001 
                # it is the 3rd table.
                tables = archive_soup.find_all("table")
                table = tables[-1]
                # The tables generally have 3 columns for - Atlantic, E. Pacific, and C. Pacific.
                # We ony require Atlantic. Find the column for Atlantic.
                i_atlantic = 0
                for col_header in table.find_all("th"):
                    if col_header.get_text() == "Atlantic":
                        break
                    else:
                        i_atlantic += 1
                # Get the list of Atlantic storms for each year.
                storm_list = table.find_all("td")[i_atlantic]
                # From the list, access the web page for each storm using the url for each storm.
                for storm_url in storm_list.find_all("a"):
                    # From the list of storms filter out hurricanes.
                    if "HURRICANE" in storm_url.get_text() or "Hurricane" in storm_url.get_text():
                        storm_name = storm_url.get_text().replace(" ", "_").lower()
                        print("[{}]".format(datetime.now().time().strftime("%H:%M")), year, storm_name)
                        storm_url = url_for_each_year+storm_url.get("href")
                        time_stamps, storm_path, pressures, wind_speeds, predicted_paths = self.__generate_storm_track(year, storm_url)
                        # Append to file           
                        f.write("{year}\t{storm_name}\t{time_stamps}\t{storm_path}\t{pressures}\t{wind_speeds}\t{predicted_paths}\n".format(
                            year = year, 
                            storm_name = storm_name, 
                            time_stamps = time_stamps, 
                            storm_path = storm_path, 
                            pressures = pressures, 
                            wind_speeds = wind_speeds, 
                            predicted_paths = predicted_paths))
            f.close()
        

    def __convert_latitude_longitude(self, latitude, longitude):
        '''
        Helper function to set the correct sign for latitudes and longitudes
        '''
        lat = 0
        long = 0
        # If latitude contains "N", then remove it and convert it to a +ve float, else remove "S" and convert it to a -ve float.
        if "N" in latitude:
            lat = float(latitude.replace("N", ""))
        else:
            lat = -1 * float(latitude.replace("S", ""))
        # If longitude contains "W", then remove it and convert it to a -ve float, else remove "E" and convert it to a +ve float.
        if "W" in longitude:
            long = -1 * float(longitude.replace("W", ""))
        else:
            long = float(longitude.replace("E", ""))
        return (lat, long)

    def __get_html(self, url):
        '''
        Helper function to get the html text from a url.
        '''
        html = None
        while html == None:
            try:
                html = requests.get(url).text
            except Exception as e:
                print("Attempt to access url {} failed with the message {}. Reattempting...".format(url, str(e)))
                continue
        return html

    def __line_to_list(self, line):
        '''
        Helper function to split a line to list of words.
        '''
        words = [word.strip() for word in line.split()]
        words = list(filter(None, words)) # Remove any empty strings
        return words

    def __generate_storm_track(self, year, storm_url):
        '''
        Helper function to gather data for a storm. 
        '''
        storm_html = self.__get_html(storm_url)
        storm_soup = BeautifulSoup(storm_html, "html.parser")
        # Find all forecast advisories. Links to forecast advisories contains the expression "MAL" or ".fstadv.".
        forecast_links = storm_soup.find_all("a", href=[re.compile("MAL"), re.compile(".fstadv.")])
        time_stamps = ""
        storm_path = ""
        predicted_paths = ""
        pressures = ""
        wind_speeds = ""
        for link in forecast_links:
            predicted_path = ""
            link_text = link.get_text()
            index_and_time_on_link = self.__line_to_list(link_text)
            time_on_link = index_and_time_on_link[1][:4]
            forecast_url = "https://www.nhc.noaa.gov{}".format(link.get("href"))
            forecast_html = self.__get_html(forecast_url)
            forecast_soup = BeautifulSoup(forecast_html, "html.parser")
            message = forecast_soup.find("pre").text.splitlines()
            line_index = 0
            for line in message:
                try:
                    line = line.strip()
                    if line.startswith(time_on_link) and line.endswith(str(year)):
                        words = self.__line_to_list(line)
                        # Get timestamp
                        if len(words) >= 5:
                            year = int(words[-1])
                            day = int(words[-2])
                            month_name = words[-3]
                            time = words[0]
                            hour = int(time[:2])
                            mins = int(time[2:4])
                            # Convert month name to number
                            month_number = int(datetime.strptime(month_name, "%b").month)
                            time = "({},{},{},{},{})".format(year, month_number, day, hour, mins)
                            if time_stamps == "":
                                time_stamps += time 
                            else:
                                time_stamps += ";"+time
                    if "REPEAT...CENTER LOCATED NEAR" in line:
                        words = self.__line_to_list(line)
                        # Get latitude and longitude
                        if len(words) >= 5:
                            latitude, longitude = self.__convert_latitude_longitude(words[3], words[4])
                            if storm_path == "":
                                storm_path += "({},{})".format(longitude, latitude)
                            else:
                                storm_path += "/({},{})".format(longitude, latitude)
                    if any(text in line for text in ["REPEAT...CENTER LOCATED NEAR", "FORECAST VALID", "OUTLOOK VALID"]):
                        # Sometimes line may contain additional information append at the end using the expression ...
                        # eg: OUTLOOK VALID 02/0600Z 47.0N  70.0W...EXTRATROPICAL
                        # eg: FORECAST VALID 30/0600Z 32.3N  89.3W...INLAND
                        # eg: OUTLOOK VALID 14/1800Z...DISSIPATED
                        # eg: OUTLOOK VALID 02/1800Z...ABSORBED
                        # These additional information are not necessary and so remove it. 
                        # However, remove care fully, since the expression ... is also associated with REPEAT..., which should not be removed.
                        if "REPEAT" not in line:
                            # Remove the additional information presented after "..."
                            line = line.split("...")[0]
                        line = line[:-1] if line.endswith(".") else line
                        words = self.__line_to_list(line)
                        # Get latitude and longitude
                        if len(words) >= 5:
                            latitude, longitude = self.__convert_latitude_longitude(words[3], words[4])
                            if predicted_path == "":
                                predicted_path += "[({},{})".format(longitude,latitude)
                            else:
                                predicted_path += "/({},{})".format(longitude,latitude)
                    if "ESTIMATED MINIMUM CENTRAL PRESSURE" in line:
                        words = self.__line_to_list(line)
                        # Get pressure
                        if len(words) >= 5:
                            pressure = words[4]
                            if pressures == "":
                                pressures += str(pressure)
                            else:
                                pressures += ";{}".format(pressure)
                    if "MAX SUSTAINED WINDS" in line:
                        words = self.__line_to_list(line)
                        # Get wind speed.
                        if len(words) >= 4:
                            wind_speed = words[3]
                            if wind_speeds == "":
                                wind_speeds += str(wind_speed)
                            else:
                                wind_speeds += ";{}".format(wind_speed)  
                    # Increment line index
                    line_index += 1
                except Exception as e:
                    print("Error {} while parsing \"{}\". Link for {} {}".format(str(e), line, link_text, link))
            predicted_path += "]"
            if predicted_paths == "":
                predicted_paths += predicted_path
            else:
                predicted_paths += ";"+predicted_path
        return (time_stamps, storm_path, pressures, wind_speeds, predicted_paths)

    def get_storms(self):
        '''
        Returns the storms archive as a Pandas dataframe.
        '''
        f = open(self.file_path, "r")
        df = pd.DataFrame(columns=["year", "storm_name", "time_stamps", "path(long,lat)", "min_atm_pressure(MB)", "max_sustained_windspeed(knots)", "predicted_paths(long,lat)"])
        lines = f.readlines()
        for i in range(1,len(lines)):
            line = lines[i]
            year, storm_name, time_stamps, path, pressures, wind_speeds, predicted_paths = line.split()
            year = float(year)
            df_time_stamps = []
            for time in time_stamps.split(";"):
                time = datetime.strptime(time, "(%Y,%m,%d,%H,%M)")
                df_time_stamps.append(time)
            df_path = []
            for location in path.split("/"):
                location = location[1:-1] # remove the brackets
                long, lat = location.split(",")
                df_path.append((float(long), float(lat)))
            df_pressures = [float(pressure) for pressure in pressures.split(";")]
            df_wind_speeds = [float(wind_speed) for wind_speed in wind_speeds.split(";")]
            df_predicted_paths = []
            for path in predicted_paths.split(";"):
                path = path[1:-1] # remove brackets
                points = []
                path = list(filter(None, path.split("/"))) # Remove any empty strings after splitting. 
                for point in path:
                    point = point[1:-1] # remove brackets
                    long, lat = point.split(",")
                    points.append((float(long), float(lat)))
                df_predicted_paths.append(points)
            # Append to df
            df.loc[len(df)] = [year, storm_name, df_time_stamps, df_path, df_pressures, df_wind_speeds, df_predicted_paths] 
        f.close()
        return df     

        
if __name__ == '__main__':   
    module_dir = pathlib.Path(__file__).parent.resolve()
    root_dir = module_dir.parent.parent
    db_file = root_dir.joinpath("data", "storms")
    storms = Storms_archieve(db_file)
    df = storms.get_storms()
    print(df)