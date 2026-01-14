You can compile `VGainScans_ana.cpp`  and then quickly run analyses in parallel using GNU Parallel with the following command:

```bash
parallel -j 10 ./launch_ana.sh ::: config/config_module_*
```
