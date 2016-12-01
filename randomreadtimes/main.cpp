#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include "picosha2.h"

using namespace std;

char get_value(int global_index) {
	int hash_id = global_index/64;
	int index = global_index%64;
	
	stringstream strs;
	strs << hash_id;
	string hashed_str;
	picosha2::hash256_hex_string(strs.str(), hashed_str);

	return hashed_str[index];
}

void flush_all_caches() {
	const int size = 100*1024*1024; // Allocate 100M. Set much larger then L3
	char *c = (char *)malloc(size);
	for (int i = 0; i < 0x3f; i++)
		for (int j = 0; j < size; j++)
			c[j] = i*j;
	free(c);
	system("sync");
	system("echo 3 | sudo /usr/bin/tee /proc/sys/vm/drop_caches");
	system("sync");
}

string mem_database;

int main() {
	int mem_database_size = 20000000;
	
	cout << "RAND_MAX..." << RAND_MAX << " should be > " << mem_database_size << endl;
	cout << "Filling in memory database..." << endl;
	stringstream mem_database_stream;
	stringstream strs;
	for(int i=0; i<mem_database_size; i++) {
		stringstream strs;
		strs << i;
		string hashed_str;
		picosha2::hash256_hex_string(strs.str(), hashed_str);
		mem_database_stream << hashed_str;
	}
	mem_database = mem_database_stream.str();
	
	cout << "Size of database: " << mem_database.size()/1000 << "KB" << "   (" << (mem_database.size()/1000)/(8192.0) << " x L3)" << endl;

	int n_exp_raw = mem_database_size/40000.0;//40K hopefully >> size of cacheline

	for(double factor = 10.0/(double)mem_database_size; factor<0.1; factor*=2) {
		int n_exp = mem_database_size*factor;
		n_exp = 10000000;
		
		ifstream fin("database.txt");
		cout << "Running " << n_exp << " reads" << endl;
		cout << "------------------" << endl;





		int serialindices[n_exp];
		serialindices[0] = rand()%(mem_database_size-n_exp-1);
		for(int i=0; i<n_exp; i++)
			serialindices[i] = serialindices[i-1]+1;


		cout << "Checking serial baseline..." << endl;

		flush_all_caches();
		auto t_serial_base_begin = chrono::high_resolution_clock::now();
		int total_serial_base = 0;
		for(int i=0; i<n_exp; i++) {
			total_serial_base += serialindices[i];
		}
		auto t_serial_base_end = chrono::high_resolution_clock::now();
		auto t_serial_base = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_serial_base_end-t_serial_base_begin).count());

		cout << "Checking serial in-memory..." << endl;
			
		flush_all_caches();
		auto t_serial_mem_begin = chrono::high_resolution_clock::now();
		int total_serial_mem = 0;
		for(int i=0; i<n_exp; i++) {
			total_serial_mem += mem_database[serialindices[i]];
		}
		auto t_serial_mem_end = chrono::high_resolution_clock::now();
		auto t_serial_mem = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_serial_mem_end-t_serial_mem_begin).count());

		//cout << "Checking serial disc..." << endl;

		//auto t_serial_disc_begin = chrono::high_resolution_clock::now();
		//int total_serial_disc = 0;
		//for(int i=0; i<n_exp; i++) {
		//	fin.seekg((streampos) serialindices[i]);
		//	total_serial_disc += (char)fin.peek();
		//}
		//auto t_serial_disc_end = chrono::high_resolution_clock::now();
		//auto t_serial_disc = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_serial_disc_end-t_serial_disc_begin).count());

		cout << "Checking serial disc..." << endl;

		flush_all_caches();
		auto t_serial_disc_begin = chrono::high_resolution_clock::now();
		int total_serial_disc = 0;
		fin.seekg((streampos) serialindices[0]);
		for(int i=0; i<n_exp; i++) {
			total_serial_disc += fin.get();
		}
		auto t_serial_disc_end = chrono::high_resolution_clock::now();
		auto t_serial_disc = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_serial_disc_end-t_serial_disc_begin).count());









		int indices[n_exp];
		for(int i=0; i<n_exp; i++)
			indices[i] = rand()%mem_database_size;

		
		cout << "Checking baseline..." << endl;

		flush_all_caches();
		auto t_base_begin = chrono::high_resolution_clock::now();
		int total_base = 0;
		for(int i=0; i<n_exp; i++) {
			total_base += indices[i];
		}
		auto t_base_end = chrono::high_resolution_clock::now();
		auto t_base = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_base_end-t_base_begin).count());

		cout << "Checking in-memory..." << endl;
			
		flush_all_caches();
		auto t_mem_begin = chrono::high_resolution_clock::now();
		int total_mem = 0;
		for(int i=0; i<n_exp; i++) {
			total_mem += mem_database[indices[i]];
		}
		auto t_mem_end = chrono::high_resolution_clock::now();
		auto t_mem = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_mem_end-t_mem_begin).count());

		cout << "Checking disc..." << endl;

		flush_all_caches();
		auto t_disc_begin = chrono::high_resolution_clock::now();
		int total_disc = 0;
		for(int i=0; i<n_exp; i++) {
			fin.seekg((streampos) indices[i]);
			total_disc += (char)fin.peek();
		}
		auto t_disc_end = chrono::high_resolution_clock::now();
		auto t_disc = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_disc_end-t_disc_begin).count());




		int do_not_optimize = total_disc + total_mem + total_base +
				  total_serial_disc + total_serial_mem + total_serial_base;

		cout << "IO of result to avoid over-optimization by gcc: " << do_not_optimize << endl;

		cout << "Results: " << endl;
		cout << "\"base\"," << t_base << "/" << n_exp << "," << n_exp << ","<< endl;
		cout << "\"mem\"," << t_mem  <<  "/" << n_exp << "," << n_exp << ","<< endl;
		cout << "\"disc\"," << t_disc << "/" << n_exp << "," << n_exp << ","<< endl;
		cout << "\"baseS\"," << t_serial_base << "/" << n_exp << "," << n_exp << "," << endl;
		cout << "\"memS\"," << t_serial_mem   << "/" << n_exp << "," << n_exp << "," << endl;
		cout << "\"discS\"," << t_serial_disc << "/" << n_exp << "," << n_exp << "," << endl;
		break;//!
	}

	return 0;
}
		//cout << "calculated: " << get_value(indices[i]) << endl;
		//cout << "in memory:  " << mem_database[indices[i]] << endl;
		//fin.seekg((streampos) indices[i]);
		//cout << "on disc:    " << (char)fin.peek() << endl;
