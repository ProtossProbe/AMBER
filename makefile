clean:
	rm assets/_output/Ast_7*.txt assets/_output/Ast_6*.txt assets/_output/*.txt assets/_output/* src/tempCodeRunnerFile.py

run:
	./main assets/input.in assets/info.in assets/switch.in

plot:
	python src/python/_analytical.py

compile:
	g++-7 -g -fopenmp -std=c++17 -Wall -Wno-unused-variable -Wno-unused-function -Wno-reorder src/crtbp.cpp src/resonance.cpp src/main.cpp -o main