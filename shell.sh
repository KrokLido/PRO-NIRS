#!/bin/bash

# This is a comment
echo "Shell scrypt run..."

#python3.11 3d_main.py to_point
#python3.11 3d_main.py dog
python3.11 3d_main.py fields
#data_from_srw='python3.11 3d_main.py to_point'
#echo data_from_srw

echo "Gnuplot..."

# График координат в пространстве.
echo "Graphic xyz..."
#gnuplot graphic_xyz.gpi
rm DUMP_graphic_xyz.txt

# График координат масс.
echo "Graphic m..."
#gnuplot graphic_m.gpi
rm graphic_m.txt

# График координат скоростей в пространстве.
echo "Vx, Vy, Vz..."
gnuplot graphic_v.gpi
gnuplot graphic_vsum.gpi
gnuplot graphic_fields.gpi
#rm graphic_vx.txt
#rm graphic_vy.txt
#rm graphic_vz.txt
#okular

echo "End of shell-code."
