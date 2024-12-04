#!/bin/bash

# Please do not modify this file. Refer to the comments presented
# in docker-mitoseg.sh (Linux/Mac) or docker-mitoseg.cmd (for
# Microsoft Windows) files for customizing the utility execution.

/app/build/mitoseg --src /app/data/src --dst /app/data/output --settings-file /app/data/settings.yaml "$@"
