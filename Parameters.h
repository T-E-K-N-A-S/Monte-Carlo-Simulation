#pragma once

//* DEFINE ALL THE PARAMETERS HERE */
#define DIM         2		// dimension of lattice // here 2; square lattice
#define V			N+2     // source + width of edge matrix + sink

// new
#define LAT 15		/*	lattice size			*/
#define N LAT+2			/*	augumented lattice size	*/
#define NDIM 2			/*	dimension of lattice	*/
#define NS pow(LAT,NDIM)	/*	*/
#define NCONFIG 1

/*	MC	*/
#define TEMP 10
#define MEAS 15 
#define SKIP 500

/*	cluster	*/
#define INFINITE 10000
#define HIGH INFINITE
#define VER N-1
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))

//delta is multipling factor with Bmat
#define delta		1
#define del_beg     .7
#define del_end     1.2
#define del_inc     .1

 
// text substitute
#define tab			"\t"
#define d			" : "

// used in push relabel
#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))