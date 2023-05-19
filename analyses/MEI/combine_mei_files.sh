mei_path=/global/scratch/users/mvazquez/projects/annocomba/results/
(for s in mMyoAui1 mMyoCai1 mMyoEvo1 mMyoLuc1 mMyoOcc1 mMyoThy1 mMyoVel1 mMyoVol1 mMyoYum1
do 
    path=./input/${mei_path}/${s}/REPEATMASKER/full/${s}_sorted.fas.out
    #echo $path
    >&2 echo $path
    tail -n +4 $path | awk '{$1=$1};{if (NF>14){print $0}}' | cut -d " " -f 1-15 | awk -v s=$s -v OFS=" " '{print s " " $0}' | tr -s " " 
done) | gzip -c >input/combined_mei.out.gz

