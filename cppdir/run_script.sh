#!/bin/bash
echo "executing simulation..."
./main.exe > out.txt
echo "simulation finished. Launching plotting program..."
python -i pyplotter.py
