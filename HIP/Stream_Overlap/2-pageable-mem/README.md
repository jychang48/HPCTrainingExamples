This Stream Overlap example uses pageable memory H2D data copy across streams.

## Build and run
```
cd /path/to/Stream_Overlap/2-pageable-mem
mkdir build
cd build
cmake ../
make -j

./compute_comm_overlap <num-of-streams>

# Unfied shared memory run
HSA_XNACK=1 ./compute_comm_overlap <num-of-streams>
```

## Profile using Omnitrace (1.11.3)
```
omnitrace-instrument -o compute_comm_overlap.inst -- compute_comm_overlap
omnitrace-run -- ./compute_comm_overlap.inst <num-of-streams>
```
