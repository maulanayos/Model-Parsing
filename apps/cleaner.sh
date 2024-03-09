#!/bin/bash
find /mnt/data/input/ecmwf_v2/unzipped/2024 -type d -mtime +2 -not -path '/mnt/data/input/ecmwf_v2/unzipped/2024' -exec rm -r {} \;
find /mnt/data/input/ecmwf_v2/zip/2024 -type d -mtime +2 -not -path '/mnt/data/input/ecmwf_v2/zip/2024' -exec rm -r {} \;
find /mnt/data/input/arpege/2024 -type d -mtime +2 -not -path '/mnt/data/input/arpege/2024' -exec rm -r {} \;
find /mnt/data/input/arome/2024 -type d -mtime +1 -not -path '/mnt/data/input/arome/2024' -exec rm -r {} \;
