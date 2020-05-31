#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vld.h>
#include "surf.h"
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
char *inputString(FILE* fp, size_t size)
{
	//The size is extended by the input with the value of the provisional
	char *str;
	int ch;
	size_t len = 0;
	str = realloc(NULL, sizeof(char)*size);//size is start size
	if (!str)return str;
	while (EOF != (ch = fgetc(fp)) && ch != '\n') {
		str[len++] = ch;
		if (len == size) {
			str = realloc(str, sizeof(char)*(size += 16));
			if (!str)return str;
		}
	}
	str[len++] = '\0';
	return realloc(str, sizeof(char)*len);
}
#include <math.h>
void setPixel(float *mat, float fx, float fy, float r, float g, float b, int width, int height, int colourComponent, int is_added)
{
	const int
		x = (int)fx - (fx >= 0 ? 0 : 1), nx = x + 1,
		y = (int)fy - (fy >= 0 ? 0 : 1), ny = y + 1;
	const float
		dx = fx - x,
		dy = fy - y;
	if (y >= 0 && y < height)
	{
		if (x >= 0 && x < width)
		{
			const float w1 = (1.0f - dx)*(1.0f - dy);
			float w2 = is_added ? 1.0f : (1.0f - w1);
			int pos = y * width * colourComponent + x * colourComponent;
			mat[pos] = (w1 * r) + (w2 * mat[pos]);
			mat[pos + 1] = (w1 * g) + (w2 * mat[pos + 1]);
			mat[pos + 2] = (w1 * b) + (w2 * mat[pos + 2]);
		}
		if (nx >= 0 && nx < width)
		{
			const float w1 = dx * (1.0f - dy);
			float w2 = is_added ? 1.0f : (1.0f - w1);
			int pos = y * width * colourComponent + nx * colourComponent;
			mat[pos] = (w1 * r) + (w2 * mat[pos]);
			mat[pos + 1] = (w1 * g) + (w2 * mat[pos + 1]);
			mat[pos + 2] = (w1 * b) + (w2 * mat[pos + 2]);
		}
	}
	if (ny >= 0 && ny < height)
	{
		if (x >= 0 && x < width)
		{
			const float w1 = (1.0f - dx)*dy;
			float w2 = is_added ? 1.0f : (1.0f - w1);
			int pos = ny * width * colourComponent + x * colourComponent;
			mat[pos] = (w1 * r) + (w2 * mat[pos]);
			mat[pos + 1] = (w1 * g) + (w2 * mat[pos + 1]);
			mat[pos + 2] = (w1 * b) + (w2 * mat[pos + 2]);
		}
		if (nx >= 0 && nx < width)
		{
			const float w1 = dx * dy;
			float w2 = is_added ? 1.0f : (1.0f - w1);
			int pos = ny * width * colourComponent + nx * colourComponent;
			mat[pos] = (w1 * r) + (w2 * mat[pos]);
			mat[pos + 1] = (w1 * g) + (w2 * mat[pos + 1]);
			mat[pos + 2] = (w1 * b) + (w2 * mat[pos + 2]);
		}
	}
}
void DrawCircle(float *img, float x0, float y0, float r, int width, int height, int colourComponent, float red, float green, float blue)
{
	float x = 0.0f, y = 0.0f;
	float r2 = r * r;
	float d2;
	x = r;
	while (y < x)
	{
		setPixel(img, x0 + x, y0 + y, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 + y, y0 + x, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 - x, y0 - y, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 - y, y0 - x, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 + x, y0 - y, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 - x, y0 + y, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 + y, y0 - x, red, green, blue, width, height, colourComponent, 1);
		setPixel(img, x0 - y, y0 + x, red, green, blue, width, height, colourComponent, 1);
		d2 = x * x + y * y;
		if (d2 > r2)
			x--;
		else
			y++;
	}
}
void DrawLine(float *img, float x0, float y0, float length, float ori, int width, int height, int colourComponent, float red, float green, float blue)
{
	float sinori = sinf(ori);
	float cosori = cosf(ori);
	int i;
	float x, y;
	float rouneD = roundf(length);
	for (i = 1; i < (int)rouneD; i++)
	{
		x = cosori * i;
		y = sinori * i;
		setPixel(img, x0 + x, y0 + y, red, green, blue, width, height, colourComponent, 1);
	}
	if (fabsf(length - rouneD) > 0.0f)
	{
		x = cosori * (length - 1.0f);
		y = sinori * (length - 1.0f);
		setPixel(img, x0 + x, y0 + y, red, green, blue, width, height, colourComponent, 1);
	}
	setPixel(img, x0, y0, 0.775f, 0.775f, 0.775f, width, height, colourComponent, 1);
}
#ifndef M_PI
#define M_PI 3.141592653589793f
#endif
#ifndef M_PI_2
#define M_PI_2 1.570796326794897f
#endif
float orientationCal(float p1X, float p1Y, float p2X, float p2Y)
{
	float y = p2Y - p1Y, x = p2X - p1X;
	float result;
	if (x != 0.0f)
	{
		const union { float flVal; unsigned int nVal; } tYSign = { y };
		const union { float flVal; unsigned int nVal; } tXSign = { x };
		if (fabsf(x) >= fabsf(y))
		{
			union { float flVal; unsigned int nVal; } tOffset = { M_PI };
			tOffset.nVal |= tYSign.nVal & 0x80000000u;
			tOffset.nVal *= tXSign.nVal >> 31;
			result = tOffset.flVal;
			const float z = y / x;
			result += (0.97239411f + -0.19194795f * z * z) * z;
		}
		else
		{
			union { float flVal; unsigned int nVal; } tOffset = { M_PI_2 };
			tOffset.nVal |= tYSign.nVal & 0x80000000u;
			result = tOffset.flVal;
			const float z = x / y;
			result -= (0.97239411f + -0.19194795f * z * z) * z;
		}
	}
	else if (y > 0.0f)
		result = M_PI_2;
	else if (y < 0.0f)
		result = -M_PI_2;
	return result;
}
float lengthCal(float p1X, float p1Y, float p2X, float p2Y)
{
	float x = (p1X - p2X) * (p1X - p2X) + (p1Y - p2Y) * (p1Y - p2Y);
	unsigned int i = *(unsigned int*)&x;
	i += 127 << 23;
	i >>= 1;
	return *(float*)&i;
}
void drawPoint(float *img, int width, int height, int colourComponent, Ipoint *ipt, float red, float green, float blue)
{
	float o = ipt->orientation;
	DrawCircle(img, ipt->x, ipt->y, 3.0f, width, height, colourComponent, red, green, blue);
}
void ShowKeyPointOri(SURFDescriptor *desp, float *dst, int width, int height, int colourComponent)
{
	float s, o, r1, c1;
	int lap;
	for (int i = 0; i < desp->hessians->interestPtsLen; i++)
	{
		SURFDescriptor *it = desp;
		s = (2.5f * it->hessians->ipts[i].scale);
		o = it->hessians->ipts[i].orientation;
		lap = it->hessians->ipts[i].laplacian;
		r1 = it->hessians->ipts[i].y;
		c1 = it->hessians->ipts[i].x;
		DrawLine(dst, c1, r1, s, o, width, height, colourComponent, 0.0f, 1.0, 0.0f);
		if (lap == 1)
			DrawCircle(dst, c1, r1, s, width, height, colourComponent, 0.0f, 0.0f, 1.0f);
		else
			DrawCircle(dst, c1, r1, s, width, height, colourComponent, 1.0f, 0.0f, 0.0f);
	}
}
// Descriptor visualization
/*int main()
{
	int i;
	float tclock;
	clock_t start, end;
	int bpp, width, height;
	unsigned char *in_img1 = stbi_load("sf.jpg", &width, &height, &bpp, 3);
	if (!in_img1)
	{
		printf("Image load failed\n");
		return -1;
	}
	unsigned char **deinterleaved = channel_split(in_img1, width, height, 3);
	Matrix mat = LoadToMatrixGreyWeighted(width, height, deinterleaved);
	for (i = 0; i < 3; i++)
		free(deinterleaved[i]);
	free(deinterleaved);
	SURFDescriptor surf;
	SURFInitialize(&surf, &mat, 5, 2, 0.0004f);
	start = clock();
	SURFCompute(&surf, mat.data);
	end = clock();
	tclock = (float)(end - start) / CLOCKS_PER_SEC;
	printf("JamesSURF found: %d interest points \n", surf.hessians->interestPtsLen);
	printf("JamesSURF took: %f seconds \n", tclock);
	printf("Stress test: second run\n");
	start = clock();
	SURFCompute(&surf, mat.data);
	end = clock();
	tclock = (float)(end - start) / CLOCKS_PER_SEC;
	printf("JamesSURF found: %d interest points \n", surf.hessians->interestPtsLen);
	printf("JamesSURF took: %f seconds \n", tclock);
	float *imageFloat = (float*)malloc(width * height * 3 * sizeof(float));
	for (i = 0; i < width * height * 3; i++)
		imageFloat[i] = in_img1[i] / 255.0f;
	ShowKeyPointOri(&surf, imageFloat, width, height, 3);
	for (i = 0; i < width * height * 3; i++)
	{
		float value = imageFloat[i] * 255.0f;
		if (value > 255.0f)
			in_img1[i] = 255U;
		else if (value < 0.0f)
			in_img1[i] = 0U;
		else
			in_img1[i] = (unsigned char)value;
	}
	free(imageFloat);
	SURFFree(&surf);
	free(mat.data);
	stbi_write_png("output.jpg", width, height, 3, in_img1, 3 * width);
	free(in_img1);
	return 0;
}*/
double RandomFloat(double a, double b, double amplitude)
{
	double random = ((double)rand()) / 32768.0;
	double diff = b - a;
	double r = random * diff;
	return amplitude * (a + r);
}
// Image matching
int main()
{
	srand(2018);
	clock_t start, end;
	int bpp, width1, height1, width2, height2;
	unsigned char *in_img1 = stbi_load("IMG_0561.jpg", &width1, &height1, &bpp, 3);
	unsigned char *in_img2 = stbi_load("IMG_0562.jpg", &width2, &height2, &bpp, 3);
	if (!in_img1 || !in_img2)
	{
		printf("Image load failed\n");
#ifdef _WIN32
		system("pause");
#endif
		return -1;
	}
	unsigned char **deinterleaved1 = channel_split(in_img1, width1, height1, 3);
	Matrix mat1 = LoadToMatrixGreyWeighted(width1, height1, deinterleaved1);
	unsigned char **deinterleaved2 = channel_split(in_img2, width2, height2, 3);
	Matrix mat2 = LoadToMatrixGreyWeighted(width2, height2, deinterleaved2);
	for (int i = 0; i < 3; i++)
	{
		free(deinterleaved1[i]);
		free(deinterleaved2[i]);
	}
	free(deinterleaved1);
	free(deinterleaved2);
	SURFDescriptor surf1, surf2;
	SURFInitialize(&surf1, &mat1, 4, 4, 0.0004f);
	SURFInitialize(&surf2, &mat2, 4, 4, 0.0004f);
	start = clock();
	SURFCompute(&surf1, mat1.data);
	SURFCompute(&surf2, mat2.data);
	IpointPair matches = getMatches(&surf1, &surf2, 0.75f);
	end = clock();
	float tclock = (float)(end - start) / CLOCKS_PER_SEC;
	printf("JamesSURF found: %d matches \n", matches.count);
	printf("JamesSURF took: %f seconds \n", tclock);
	SURFFree(&surf1);
	SURFFree(&surf2);
	free(mat1.data);
	free(mat2.data);
	float *imageFloat1 = (float*)malloc(width1 * height1 * 3 * sizeof(float));
	float *imageFloat2 = (float*)malloc(width2 * height2 * 3 * sizeof(float));
	int i;
	for (i = 0; i < width1 * height1 * 3; i++)
		imageFloat1[i] = in_img1[i] / 255.0f;
	for (i = 0; i < width2 * height2 * 3; i++)
		imageFloat2[i] = in_img2[i] / 255.0f;
	for (i = 0; i < matches.count; ++i)
	{
		float red = (float)RandomFloat(0.0, 1.0, 1.0);
		float green = (float)RandomFloat(0.0, 1.0, 1.0);
		float blue = (float)RandomFloat(0.0, 1.0, 1.0);
		drawPoint(imageFloat1, width1, height1, 3, &matches.pt1[i], red, green, blue);
		drawPoint(imageFloat2, width2, height2, 3, &matches.pt2[i], red, green, blue);
		float circleOrientation = orientationCal(matches.pt1[i].x, matches.pt1[i].y, matches.pt2[i].x + width1, matches.pt2[i].y);
		float distance = lengthCal(matches.pt1[i].x, matches.pt1[i].y, matches.pt2[i].x + width1, matches.pt2[i].y);
		DrawLine(imageFloat1, matches.pt1[i].x, matches.pt1[i].y, distance, circleOrientation, width1, height1, 3, red, green, blue);
		circleOrientation = orientationCal(matches.pt1[i].x - width1, matches.pt1[i].y, matches.pt2[i].x, matches.pt2[i].y);
		distance = lengthCal(matches.pt1[i].x - width1, matches.pt1[i].y, matches.pt2[i].x, matches.pt2[i].y);
		DrawLine(imageFloat2, matches.pt1[i].x - width1, matches.pt1[i].y, distance, circleOrientation, width2, height2, 3, red, green, blue);
	}
	matcherFree(&matches);
	for (i = 0; i < width1 * height1 * 3; i++)
	{
		float value = imageFloat1[i] * 255.0f;
		if (value > 255.0f)
			in_img1[i] = 255U;
		else if (value < 0.0f)
			in_img1[i] = 0U;
		else
			in_img1[i] = (unsigned char)value;
	}
	free(imageFloat1);
	for (i = 0; i < width2 * height2 * 3; i++)
	{
		float value = imageFloat2[i] * 255.0f;
		if (value > 255.0f)
			in_img2[i] = 255U;
		else if (value < 0.0f)
			in_img2[i] = 0U;
		else
			in_img2[i] = (unsigned char)value;
	}
	free(imageFloat2);
	stbi_write_png("output1.png", width1, height1, 3, in_img1, 3 * width1);
	stbi_write_png("output2.png", width2, height2, 3, in_img2, 3 * width2);
	free(in_img1);
	free(in_img2);
	return 0;
}