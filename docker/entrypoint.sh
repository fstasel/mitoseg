#!/bin/bash

/app/build/mitoseg -src /app/data/src -dst /app/data/results -settingsFile /app/data/settings.yaml "$@"
