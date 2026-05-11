all:
	g++ -O3 -march=native -flto=auto -pthread fast39.cc -o fast39
	g++ -O3 -o makesievebase makesievebase.cc intpoly.cc mpz_poly.cc factorsmall.cc -std=c++11 -fopenmp -lgmp -lgmpxx
	g++ -O3 -o hdlcs hdlcs.cc mpz_poly.cc factorsmall.cc intpoly.cc -lquadmath -fext-numeric-literals -std=c++11 -pthread -lgmp -lgmpxx
	#g++ -O3 -march=native -flto=auto -fno-fat-lto-objects -fuse-linker-plugin -ffp-contract=fast fast38.cc -o fast38

debug:
	g++ -O0 -g fast39.cc -o fast39
	g++ -O0 -g -o makesievebase makesievebase.cc intpoly.cc mpz_poly.cc factorsmall.cc -std=c++11 -fopenmp -lgmp -lgmpxx
	g++ -O0 -g -o hdlcs hdlcs.cc mpz_poly.cc factorsmall.cc intpoly.cc -lquadmath -fext-numeric-literals -std=c++11 -pthread -lgmp -lgmpxx -fsanitize=address -fsanitize=undefined

clean:
	rm fast39
	rm makesievebase
	rm hdlcs
