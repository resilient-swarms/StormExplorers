This is the standard bundle for National Oceanographic Data Center (NODC)
accession data.  Each accession maintained by NODC is given a unique integer
identifier known as an 'accession id'.  All information related to that
accession is contained within a directory using the accession id as the
directory name.  The accession directory has the following standard structure:

<accession id>/<version>-version: This directory.  Contains the <version>th
         version of the information for this accession.  The first version is
         "01-version".  Any time a new version of the originator's data is
         received, a new <version>-version directory is created.

  NODC-Readme.txt: This file.

  about: Directory.  Contains all accession related metadata including but
         not limited to the following standard file:

    journal.txt: Text file.  Contains any notes, correspondence, etc.,
         relating to this accession.

  data: Directory.  All accession data is located in the 'data' directory.

    0-data: Directory.  Contains the originator's data unmodified from its
         initial digital format as submitted to NODC.  The initial source for
         this data should be documented in the header of the
         <accession id>/<version>-version/about/journal.txt file after the
         'Source' keyword.

    1-data: Optional directory.  May contain processed version of originator's
         data from '0-data' directory, e.g., unzipped, uncompressed, untarred,
         or otherwise extracted or modified data.  A note should be found in
         <accession id>/<version>-version/about/journal.txt explaining how
         files in 1-data were derived from the files in 0-data.

    <n>-data: Optional directories.  Additional processed forms of originators
         data.  Similar to 1-data above.

<accession id>/<accession id>.<version>-version.md5: Text file.  Contains MD5
         checksums for all files in this version of the accession (except for
         that of the <accession id>.<version>-version.md5 file itself).


For further information about this accession see:

./about/journal.txt
