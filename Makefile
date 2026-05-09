all:
	g++ -O3 -march=native -flto=auto -pthread fast39.cc -o fast39
	#g++ -O3 -march=native -flto=auto -fno-fat-lto-objects -fuse-linker-plugin -ffp-contract=fast fast38.cc -o fast38

debug:
	g++ -O0 -g fast39.cc -o fast39

clean:
	rm fast39
