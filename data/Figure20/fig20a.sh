cd ../../build/Release

chmod +x isosurfacing

./isosurfacing ../../data/Figure20/grid_1.json ../../data/Figure20/csg_examples_1.json -t 0.01 -o "CSG" --tree ../../data/Figure20/csg_examples_1_tree.json
