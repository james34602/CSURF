#include "matrix.h"
typedef struct
{
	int width;
	int	height;
	int	step;
	int	filter;
	float *responses;
	unsigned char *laplacian;
} ResponseLayer;
typedef struct
{
	//! Coordinates of the detected interest point
	float x, y;
	//! Detected scale
	float scale;
	//! Orientation measured anti-clockwise from +ve x-axis
	float orientation;
	//! Sign of laplacian for fast matching purposes
	unsigned char laplacian;
	//! Array of descriptor components
	float descriptor[64]; //float descriptor[128];
	float dx, dy;
	//! Used to store cluster index
	int clusterIndex;
} Ipoint;
typedef struct
{
	Matrix *img;
	int i_width;
	int i_height;
	Ipoint *ipts;
	int iptssize;
	ResponseLayer responseMap[12];
	int octaves;
	int init_sample;
	float thresh;
	int interestPtsLen;
} FastHessian;
typedef struct
{
    Matrix *integralImg;
	FastHessian *hessians;
} SURFDescriptor;
typedef struct
{
	Ipoint *pt1, *pt2;
	int size;
	int count;
} IpointPair;
void SURFInitialize(SURFDescriptor *surf, Matrix *img, int octaves, int init_sample, float thres);
void SURFCompute(SURFDescriptor *surf, float *img);
void SURFFree(SURFDescriptor *surf);
IpointPair getMatches(SURFDescriptor *surf1, SURFDescriptor *surf2, float threshold);
void matcherFree(IpointPair *pair);