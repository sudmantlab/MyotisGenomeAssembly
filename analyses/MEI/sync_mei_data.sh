
mei_path=/global/scratch/users/mvazquez/projects/annocomba/results/
(for s in mMyoAui1  mMyoCai1  mMyoEvo1  mMyoLuc1  mMyoOcc1  mMyoThy1  mMyoVel1  mMyoVol1  mMyoYum1
do 
    path=${mei_path}/${s}/REPEATMASKER/full/${s}_sorted.fas.out
    echo $path
done)>sync_files.txt
echo rsync -arvz --progress --files-from sync_files.txt psudmant@dtn.brc.berkeley.edu:/ ./input
