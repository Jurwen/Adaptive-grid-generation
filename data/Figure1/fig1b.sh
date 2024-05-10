cd ../../build/Release

chmod +x isosurfacing

./isosurfacing ../../data/mesh/cube6.msh ../../data/Figure1/10-wikiBall.json -t 0.0005 "-o" "CSG" "--tree" ../../data/Figure1/10-wikiBall-tree.json
