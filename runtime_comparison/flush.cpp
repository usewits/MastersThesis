#include <iostream>
#include <stdlib.h>
#include <time.h>
#include <string>
#include <string.h>

using namespace std;

void flush_all_caches() {
    cout << "flushing.." << endl;

    const long long size = 270465449984LL/32LL;//Allocate 1/32 of total memory
    const long long alloc_size = 1024*1024+64;//This approximately divises size
    const int n_blocks = size/alloc_size;
    cout << "Allocate: " << size << "B = " << size/(1000.*1000.*1000.) << "G" << endl;
    cout << "In blocks:" << alloc_size << "B = " << alloc_size/(1000.*1000.) << "M" << endl;
    cout << "N blocks:" << n_blocks << endl;

    char **c = (char**)malloc(n_blocks*sizeof(char*));

    for(int k=0; k<n_blocks; k++) {
        c[k] = (char*) calloc(alloc_size, sizeof(char));
        for(int j=0; j<alloc_size; j+=64)
            c[k][j] = k*j;
    }

    unsigned int block_internal = rand();
    unsigned int index_internal = rand();
    unsigned int index2_internal = rand();
    unsigned const int flush_block_size = 1024;
    int block, index, index2;
    int result = 0;
    int n_loops = 5;
    for(int k=0; k<n_loops; k++) {
        cout << "flushing " << k << "/" << n_loops << ".." << endl;
        for(long long z=0; z<size/flush_block_size; z++) {//every byte once on avg
            block_internal = (block_internal*1103515245U+12345U);
            index_internal = (index_internal*1103515245U+12345U);
            block = (block_internal)%n_blocks;
            index = (index_internal)%(alloc_size-flush_block_size);

            int strafe = (block_internal+index_internal)%16+16;

            for(int i=0; i<flush_block_size; i+=strafe)
                c[block][index+i] += block_internal*index_internal;

            index2_internal = (index2_internal*1103515245U+12345U);
            index2 = (index2_internal)%n_blocks;
            result += c[block][index2];
        }
    }
    cout << "avoid opt. " << result << endl;
    for(int k=0; k<n_blocks; k++) {
        free(c[k]);
    }
    free(c);
    system("sync");
    system("echo 3 | sudo /usr/bin/tee /proc/sys/vm/drop_caches");
    system("sync");
    cout << "done!" << endl;
}


int main(int argc, char** argv) {
    if(argc != 2) {
        cout << "ERROR: give 1 arg" << endl;
        return 0;
    }
    int seed = 0;
    string arg(argv[1]);
    int pow = 1;
    for(int i=arg.size()-1; i>=0; i--) {
        seed += pow * (arg[i]-'0');
        pow*=10;
    }
    srand((seed+982587397)*time(NULL));
    flush_all_caches();
    cout << "!!!" << seed << " is done!!!" << endl;
    return 0;
}

