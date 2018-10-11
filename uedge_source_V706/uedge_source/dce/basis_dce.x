const MAXARRAY = 100000;
const MAXVERTS = 400000;

struct dtm_server {
		string name<256>;
};

struct image {
    int dtmport;
    string title<256>;
    float imagedata<MAXARRAY>;
    int xlen;
    int ylen;
};
struct rpcsignal {
    string params<256>;
    float xdata<MAXARRAY>;
    float signaldata<MAXARRAY>;
    int xlen;
};
struct pr_image {

    int dtmport;
    string title<256>;
    float prdata<MAXARRAY>;
    float r_verts<MAXVERTS>;
    float z_verts<MAXVERTS>;
    int xlen;
    int ylen;
    int outxlen;
    int outylen;
    float rmin;
    float rmax;
    float zmin;
    float zmax;

};
struct pr_image_int {

    int dtmport;
    string title<256>;
    float prdata<MAXARRAY>;
    float r_verts<MAXVERTS>;
    float z_verts<MAXVERTS>;
		float r_int<MAXVERTS>;
		float z_int<MAXVERTS>;
		float v_int<MAXVERTS>;
    int xlen;
    int ylen;
    int outxlen;
    int outylen;
    float rmin;
    float rmax;
    float zmin;
    float zmax;

};

typedef dtm_server dtmserver;
typedef struct image ledgeimage;
typedef struct rpcsignal ledgesignal;
typedef struct pr_image ledge_pr_image;
typedef struct pr_image_int ledge_pr_intimage;
typedef float rawimage<MAXVERTS>;


/*
The rzxform routines are now incorporated directly into
uedge instead of being an rpc.
*/
program DCEPROG {
    version DCEVERS {
        int INIT_DTM(dtmserver) = 1;
        int DTMOUT(ledgeimage) = 2;
        int DTMRZOUT(ledge_pr_image) = 3;
        int CLOSE_DTM(int) = 4;
        int XGRAPHOUT(ledgesignal) = 5;
        /*rawimage RZXFORM(ledge_pr_image) = 6;*/
        int DTMRZOUT_INT(ledge_pr_intimage) = 7;
        /*rawimage RZXFORM_INT(ledge_pr_intimage) = 8;*/
    }
    = 1;
}
= 200001;
