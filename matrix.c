#include <math.h>
#include <stdlib.h>
#include "matrix.h"
Matrix MatrixAllocate(int w, int h)
{
	Matrix mat;
	mat.width = w;
	mat.height = h;
	if (w > 0 && h > 0)
		mat.data = (float*)malloc(w*h*sizeof(float));
	else
		mat.data = 0;
	return mat;
}
Matrix* MatrixAllocatePtr(int w, int h)
{
	Matrix *mat = (Matrix*)malloc(sizeof(Matrix));
	mat->width = w;
	mat->height = h;
	if (w > 0 && h > 0)
		mat->data = (float*)malloc(w*h * sizeof(float));
	else
		mat->data = 0;
	return mat;
}
float get(Matrix *mat, int x, int y)
{
	return mat->data[y * mat->width + x];
}

inline void set(Matrix *mat, int x, int y, float v)
{
	//if (x < 0 || y < 0 || x >= width || y >= height)
	//	return;
	//else
	mat->data[y * mat->width + x] = v;
}

Matrix LoadToMatrix(int width, int height, unsigned char *to)
{
	Matrix mat = MatrixAllocate(width, height);
	for (int i = 0; i < mat.width * mat.height; i++)
		mat.data[i] = (float)to[i] / 255.0f;
	return mat;
}
Matrix LoadToMatrixGreyWeighted(int width, int height, unsigned char **buffer)
{
	Matrix mat = MatrixAllocate(width, height);
	int i, num_frames = width * height;
	float tofloat = 1.0f / 255.0f;
	for (i = 0; i < num_frames; i++)
		mat.data[i] = ((float)buffer[0][i] * 0.2989f + (float)buffer[1][i] * 0.587f + (float)buffer[2][i] * 0.114f) * tofloat;
	return mat;
}
unsigned char** channel_split(unsigned char* buffer, int width, int height, int num_channels)
{
	int i, num_frames = width * height;
	unsigned char **chan_buffers = (unsigned char**)malloc(num_channels * sizeof(unsigned char*));
	for (i = 0; i < num_channels; i++)
		chan_buffers[i] = (unsigned char*)malloc(num_frames * sizeof(unsigned char));
	int samples = num_frames * num_channels;
	for (i = 0; i < samples; i++)
		chan_buffers[(i % num_channels)][i / num_channels] = buffer[i];
	return chan_buffers;
}
void channel_join(unsigned char** chan_buffers, int num_channels, unsigned char* buffer, int num_frames)
{
	for (int i = 0; i < num_frames * num_channels; i++)
		buffer[i] = chan_buffers[i % num_channels][i / num_channels];
}
void gaussianKernel(double **gaussian, const int N, double sigma)
{
	int i, j;
	double sum = 0.0;
	for (i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			gaussian[i][j] = exp(-((i - N / 2.0)*(i - N / 2.0) + (j - N / 2.0)*(j - N / 2.0)) / (2.0*sigma*sigma));
			sum += gaussian[i][j];
		}
	}
}