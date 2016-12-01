#include <iostream>
#include <fstream>
#include <iterator>
#include <set>
#include <map>
#include <math.h>
#include <string>
#include <chrono>
#include <sys/stat.h>
#include <sys/mman.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>

#include "picosha2.h"

using namespace std;


int R1A_size = 200000000;
int R1B_size = 200000000;
int R2A_size = 2000;
int R2C_size = 2000;

unsigned long R1A_offset = 0;
unsigned long R1B_offset = R1A_offset+R1A_size;const char* R1A_mem;
unsigned long R2A_offset = R1B_offset+R1B_size;const char* R1B_mem;
unsigned long R2C_offset = R2A_offset+R2A_size;const char* R2A_mem;

const char* R2C_mem;

unsigned long filelen;
char* data_disk;

const char* R1A_disk;
const char* R1B_disk;
const char* R2A_disk;
const char* R2C_disk;

//The function mmapopen was obtained from https://github.com/hannesmuehleisen/leaderboard-demo with permission
void* mmapopen(char* filename, unsigned long &filelen, bool write) {
	int fd;
	struct stat sbuf;
	void *mapaddr;
	int fopenmode = write ? O_RDWR : O_RDONLY;
	int mmapmode = write ? PROT_WRITE : PROT_READ;

	if ((fd = open(filename, fopenmode)) == -1) {
		fprintf(stderr, "failed to open %s\n", filename);
		exit(1);
    }

    if (stat(filename, &sbuf) == -1) {
        fprintf(stderr, "failed to stat %s\n", filename);
        exit(1);
    }
    
    mapaddr = mmap(0, sbuf.st_size, mmapmode, MAP_SHARED, fd, 0);
    if (mapaddr == MAP_FAILED) {
        fprintf(stderr, "failed to mmap %s\n", filename);
        exit(1);
    }
    filelen = sbuf.st_size;
    return mapaddr;
}

//We did implement mmapclose ourselves
void mmapclose(void* mapaddr, unsigned long filelen) {
    if (munmap(mapaddr, filelen) != 0) {
        fprintf(stderr, "failed to munmap\n");
        exit(1);
    }
    return;
}




void flush_all_caches(bool close_mmap) {
	//system("./run_flush");
	if(close_mmap) {
		cout << "closing memory maps.." << endl;
		mmapclose(data_disk, filelen);
		system("sync");
	} else {
		cout << "skipped closing memory maps.." << endl;
	}
	

	cout << "flusing page cache.." << endl;
	system("echo 3 | sudo /usr/bin/tee /proc/sys/vm/drop_caches");

	system("sync");

	cout << "(re)opening memory maps.." << endl;
	filelen=-1;
	data_disk = (char*) mmapopen("database.txt", filelen, false);
	if(filelen < R1A_size+R1B_size+R2A_size+R2C_size) {
		cout << "WARNING: database is too small!" << endl;
	}
	
	R1A_disk = data_disk + R1A_offset;
	R1B_disk = data_disk + R1B_offset;
	R2A_disk = data_disk + R2A_offset;
	R2C_disk = data_disk + R2C_offset;

	system("sync");

	cout << "flusing page cache.." << endl;
	system("echo 3 | sudo /usr/bin/tee /proc/sys/vm/drop_caches");

	system("sync");


	cout << "heating up R1 and R2" << flush;
	int no_opt = 0;
	for(int i=0; i<10; i++) {
		int strafe = (rand()%16)+16;
		int woggle = (rand()%32)+32;
		for(int j=0; j<R1A_size; j+=strafe) {
			no_opt += R1A_mem[j];
			if(j%woggle == 0) {
				no_opt *= rand()+1;	
			}
		}
		for(int j=0; j<R1B_size; j+=strafe) {
			no_opt += R1B_mem[j];
			if(j%woggle == 0) {
				no_opt *= rand()+1;	
			}
		}
		cout << "." << flush;
		for(int j=0; j<R2A_size; j+=strafe) {
			no_opt += R2A_mem[j];
			if(j%woggle == 0) {
				no_opt *= rand()+1;	
			}
		}
		for(int j=0; j<R2C_size; j+=strafe) {
			no_opt += R2C_mem[j];
			if(j%woggle == 0) {
				no_opt *= rand()+1;	
			}
		}
	}
	cout << no_opt;



	system("sync");

	cout << "flusing CPU cache" << endl;

	const int size = 50*1024*1024; // Allocate 100M. Set much larger then L3
	char *c = (char *)malloc(size);
	for (int i = 0; i < 0x3f; i++)
		for (int j = 0; j < size; j++)
			c[j] = i*j;
	free(c);

	system("sync");

	cout << "wait" << endl;
	system("sleep 1");

	cout << "flusing done!" << endl;
}

double get_weight(double A, double B, double C) {
	return A+B*C;
}



char* wr_uniform_sample(const char* data, int n, int m) {//with replacement
	char* result = (char*)malloc(m*sizeof(char));
	for(int i=0; i<m; i++) {
		result[i] = data[rand()%n];
	}
	return result;
}

set<pair<int,char> > wor_uniform_sample(const char* data, int n, int m) {//without replacement
	set<pair<int,char> > result;
	while(result.size()<m) {
		int rand_index = rand()%n;
		result.insert(make_pair(rand_index, data[rand_index]));
	}
	return result;
}

//Uniform reservoir sample
char* wor_reservoir_sample(const char* data, int n, int m) {//without replacement
	char* result = (char*)malloc(m*sizeof(char));
	for(int i=0; i<m; i++) {
		result[i] = data[i];
	}
	for(int i=m; i<n; i++) {
		if(rand()%i < m) {//P(rand()%i < m) = P(U_i < m/i) = m/i, where U_i is uniform in [0,1[
			result[rand()%m] = data[i];//Replace random element
		}
	}
	return result;
}

//Based on Alg-A by Efraimidis and Spirakis
multimap<double,char> weighted_wor_reservoir_sample(const char* data, const char* w, int n, int m) {//without replacement
	double* keys = (double*)malloc(n*sizeof(double));
	for(int i=0; i<n; i++) {
		keys[i] = pow(rand()/(double)RAND_MAX,1.0/(double)w[i]);
	}

	multimap<double,char> result;
	for(int i=0; i<m; i++) {
		result.insert(make_pair(keys[i], data[i]));
	}

	for(int i=m; i<n; i++) {
		if(keys[i] > result.begin()->first) {
			auto it = result.begin();
			result.erase(it);
			result.insert(make_pair(keys[i], data[i]));
		}
	}
	free(keys);
	return result;
}

//Based on Alg-A-exp by Efraimidis and Spirakis
multimap<double,char> weighted_wor_reservoir_sample_exp(const char* data, const char* w, int n, int m) {//without replacement
	double* keys = (double*)malloc(m*sizeof(double));
	for(int i=0; i<m; i++) {
		keys[i] = pow(rand()/(double)RAND_MAX,1/(double)w[i]);
	}

	multimap<double,char> result;
	for(int i=0; i<m; i++) {
		result.insert(make_pair(keys[i], data[i]));
	}

	int i=m;
	while(true) {
		double r = rand()/(double)RAND_MAX;
		double xw = log(r)/log(result.begin()->first);
		while(xw > 0 && i < n) {
			xw -= w[i];
			i++;
		}
		if(i >= n) break;
		//At this point, xw - (w[c]+w[c+1] + ... + w[i]) <= 0 
		double tw = pow(result.begin()->first, (double)w[i]);
		double r2 = (rand()/(double)RAND_MAX)*(1-tw)+tw;
		double key = pow(r2, 1/(double)w[i]);

		auto it = result.begin();
		result.erase(it);
		result.insert(make_pair(key, data[i]));
	}
	free(keys);
	return result;
}




int main() {
  	cout << "RAND_MAX..." << RAND_MAX << " should be > " << R1A_size << endl;

	cout << "Filling in memory columns..." << endl;

	stringstream mem_database_stream;
	for(int i=0; i < (int)(R1A_size+R1B_size+R2A_size+R2C_size)/(double)8+1; i++) {
		stringstream strs;
		strs << i;
		string hashed_str;
		picosha2::hash256_hex_string(strs.str(), hashed_str);
		mem_database_stream << hashed_str;
	}
	string mem_database = mem_database_stream.str();

	R1A_mem = mem_database.c_str()+R1A_offset;
	R1B_mem = mem_database.c_str()+R1B_offset;
	R2A_mem = mem_database.c_str()+R2A_offset;
	R2C_mem = mem_database.c_str()+R2C_offset;

	cout << "Size of R1A: " << R1A_size/1000 << "KB" 
		 << "   (" << (R1A_size/1000)/(8192.0) << " x L3)" << endl;


	flush_all_caches(false);//do not close mmap (since it is not open yet)

	//Take S uniform/weighted/AWS sample in R1
	//For each element in S, take uniform/weighted sample in (part of) R2

	int do_not_optimize = 0;

	cout<<"@R1_mem"<<","<<"R2_mem"<<","<<"m"<<","<<"n1"<<","<<"n2"<<","<<"WS"<<","<<"t"<<"},{"<<endl;

	for(int experiment = 0; experiment < 4; experiment++)
	for(double m_frac = 2e-8; m_frac < 1.0; m_frac *= 10)
	for(int repeati = 0; repeati<5; repeati++)
	{
		int m = m_frac * R1A_size;
		if(m == 0) continue;
		const char* R1A;
		const char* R1B;
		const char* R2A;
		const char* R2C;
		cout << endl;
		cout << endl;
		bool R1_mem, R2_mem;
		if(experiment%2 == 0) {//0,2
			cout << "R1 on disk" << endl;
			R1A = R1A_disk;
			R1B = R1B_disk;
			R1_mem=false;
		} else {//1,3
			cout << "R1 in memory" << endl;
			R1A = R1A_mem;
			R1B = R1B_mem;
			R1_mem=true;
		}
		if(experiment/2 == 0) {//0,1
			cout << "R2 on disk" << endl;
			R2A = R2A_disk;
			R2C = R2C_disk;
			R2_mem=false;
		} else {//2,3
			cout << "R2 in memory" << endl;
			R2A = R2A_mem;
			R2C = R2C_mem;
			R2_mem=true;
		}
		cout << "m = " << m << endl;

	//WS-join, h wo c

		flush_all_caches(true);
		auto t_ws_h_wo_c_begin = chrono::high_resolution_clock::now();
		{
			multimap<double,char> S1 = weighted_wor_reservoir_sample_exp(R1A, R1B, R1A_size, m);
			vector< pair<char, char> > join_result(m);
			int index = 0;
			for(auto it : S1) {
				char* S2 = wr_uniform_sample(R2A, R2C_size, 1);
				join_result[index] = make_pair(it.second, *S2);
				do_not_optimize += (int)join_result[index].first*(int)join_result[index].second;
				index++;
				free(S2);
			}
		}
		//join_result is a sample of the join result
		auto t_ws_h_wo_c_end = chrono::high_resolution_clock::now();
		auto t_ws_h_wo_c = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_ws_h_wo_c_end-t_ws_h_wo_c_begin).count());

		cout << "WS (h w/o c) " << t_ws_h_wo_c << endl;
		cout<<"@"<<R1_mem<<","<<R2_mem<<","<<m<<","<<R1A_size<<","<<R2A_size<<","<<true<<","<<t_ws_h_wo_c<<"},{"<<endl;



	//WS-join, h wo c NO EXP

		flush_all_caches(true);
		auto t_ws_noexp_h_wo_c_begin = chrono::high_resolution_clock::now();
		{
			multimap<double,char> S1 = weighted_wor_reservoir_sample(R1A, R1B, R1A_size, m);
			vector< pair<char, char> > join_result(m);
			int index = 0;
			for(auto it : S1) {
				char* S2 = wr_uniform_sample(R2A, R2C_size, 1);
				join_result[index] = make_pair(it.second, *S2);
				do_not_optimize += (int)join_result[index].first*(int)join_result[index].second;
				index++;
				free(S2);
			}
		}
		//join_result is a sample of the join result
		auto t_ws_noexp_h_wo_c_end = chrono::high_resolution_clock::now();
		auto t_ws_noexp_h_wo_c = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_ws_noexp_h_wo_c_end-t_ws_noexp_h_wo_c_begin).count());

		cout << "WS (h w/o c) " << t_ws_noexp_h_wo_c << endl;
		cout<<"@"<<R1_mem<<","<<R2_mem<<","<<m<<","<<R1A_size<<","<<R2A_size<<","<<2<<","<<t_ws_noexp_h_wo_c<<"},{"<<endl;



	//WS-join, h w c
	//NOTE: actually O(n*m) without CDF!
//
//		flush_all_caches(true);
//		auto t_ws_h_w_c_begin = chrono::high_resolution_clock::now();
//		{
//			multimap<double,char> S1 = weighted_wor_reservoir_sample_exp(R1A, R1B, R1A_size, m);
//			vector< pair<char, char> > join_result(m);
//			int index = 0;
//			for(auto it : S1) {
//				multimap<double,char> S2 = weighted_wor_reservoir_sample_exp(R2A, R2C, R2C_size, 1);
//				join_result[index] = make_pair(it.second, S2.begin()->second);
//				do_not_optimize += (int)join_result[index].first*(int)join_result[index].second;
//				index++;
//			}
//		}
//		//join_result is a sample of the join result
//		auto t_ws_h_w_c_end = chrono::high_resolution_clock::now();
//		auto t_ws_h_w_c = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_ws_h_w_c_end-t_ws_h_w_c_begin).count());
//
//		cout << "WS (h w/ c)  " << t_ws_h_w_c << endl;



	//US-join

		flush_all_caches(true);
		auto t_us_begin = chrono::high_resolution_clock::now();
		{
			set< pair<int, char> > S1 = wor_uniform_sample(R1A, R1A_size, m);
			vector< pair<char, char> > join_result(m);
			int index = 0;
			for(auto it : S1) {
				char* S2 = wr_uniform_sample(R2A, R2C_size, 1);
				join_result[index] = make_pair(it.second, *S2);
				do_not_optimize += (int)join_result[index].first*(int)join_result[index].second;
				index++;
				free(S2);
			}
		}
		//join_result is a sample of the join result
		auto t_us_end = chrono::high_resolution_clock::now();
		auto t_us = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_us_end-t_us_begin).count());

		cout << "US           "<< t_us << endl;
		cout<<"@"<<R1_mem<<","<<R2_mem<<","<<m<<","<<R1A_size<<","<<R2A_size<<","<<false<<","<<t_us<<"},{"<<endl;


	//HWS-join

		flush_all_caches(true);
		if(m*m < R1A_size) {
			auto t_hws_begin = chrono::high_resolution_clock::now();
			{
					set< pair<int, char> > U1 = wor_uniform_sample(R1A, R1A_size, m*m);
					stringstream U1strs;
					for(auto pic : U1)
							U1strs << pic.second;
					string U1str = U1strs.str();
					const char* U1char = U1str.c_str();

					multimap<double,char> S1 = weighted_wor_reservoir_sample_exp(U1char, R1B, U1str.length(), m);
					vector< pair<char, char> > join_result(m);
					int index = 0;
					for(auto it : S1) {
							char* S2 = wr_uniform_sample(R2A, R2C_size, 1);
							join_result[index] = make_pair(it.second, *S2);
							do_not_optimize += (int)join_result[index].first*(int)join_result[index].second;
							index++;
							free(S2);
					}               
					
			}
			//join_result is a sample of the join result
			auto t_hws_end = chrono::high_resolution_clock::now();
			auto t_hws = (std::chrono::duration_cast<std::chrono::nanoseconds>(t_hws_end-t_hws_begin).count());

			cout << "HWS           "<< t_hws << endl;
			cout<<"@"<<R1_mem<<","<<R2_mem<<","<<m<<","<<R1A_size<<","<<R2A_size<<","<<3<<","<<t_hws<<"},{"<<endl;
		}


	}

	cout << "Avoid optimization: " << do_not_optimize << endl;
}

