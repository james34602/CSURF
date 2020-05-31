#ifndef _SURF_MATRIX_H_
#define _SURF_MATRIX_H_
typedef struct
{
	int width;
	int height;
	float *data;
} Matrix;
Matrix MatrixAllocate(int w, int h);
Matrix* MatrixAllocatePtr(int w, int h);
float get(Matrix *mat, int x, int y);
Matrix LoadToMatrix(int width, int height, unsigned char *to);
Matrix LoadToMatrixGreyWeighted(int width, int height, unsigned char **buffer);
unsigned char** channel_split(unsigned char* buffer, int width, int height, int num_channels);
void channel_join(unsigned char** chan_buffers, int num_channels, unsigned char* buffer, int num_frames);
void gaussianKernel(double **gaussian, const int N, double sigma);
#endif