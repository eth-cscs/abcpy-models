#!/bin/sh
cd "`dirname $0`"

# genarate BMP images
init_mono="initial.dat"
final_mono="final.dat"

init_s1="initial_S1.dat"
init_s2="initial_S2.dat"
final_s1="final_S1.dat"
final_s2="final_S2.dat"

out_init_hpc="initial_hpc.dat"
out_final_hpc="final_hpc_hpc"

head -n 2 $init_mono >  $out_init_hpc
tmp1="tmp1.txt"
tmp2="tmp2.txt"
tail -n +3 $init_s1 >$tmp1
tail -n +3 $init_s2 >$tmp2
paste $tmp1 $tmp2  -d '' >>  $out_init_hpc

head -n 2 $final_mono >  $out_final_hpc
tmp1="tmp1.txt"
tmp2="tmp2.txt"
tail -n +3 $final_s1 >$tmp1
tail -n +3 $final_s2 >$tmp2
paste $tmp1 $tmp2  -d '' >>  $out_final_hpc

./grid_to_bmp $out_init_hpc initial_hpc.bmp
./grid_to_bmp $out_final_hpc final_hpc.bmp

rm -rf "$tmp1" "$tmp2"

cd -
