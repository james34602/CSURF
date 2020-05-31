// History:
// Creation date: 2018 April
// Compute basic SURF descriptor and image matching function with multithreading
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <string.h>
#include "cpthread.h"
#include "surf.h"
int fRound(float flt)
{
	return (int)floorf(flt + 0.5f);
}
void ResLayer(ResponseLayer *Res, int width, int height, int step, int filter)
{
	Res->width = width;
	Res->height = height;
	Res->step = step;
	Res->filter = filter;
	Res->responses = (float*)malloc(width * height * sizeof(float));
	Res->laplacian = (unsigned char*)malloc(width * height * sizeof(unsigned char));
}
float get_Response(unsigned int row, unsigned int column, ResponseLayer *R)
{
	return R->responses[row * R->width + column];
}
float getResponse(unsigned int row, unsigned int column, ResponseLayer *src, ResponseLayer *R)
{
	int scale = R->width / src->width;
	return R->responses[(scale * row) * R->width + (scale * column)];
}
unsigned char getLaplacian(unsigned int row, unsigned int column, ResponseLayer *src, ResponseLayer *R)
{
	int scale = R->width / src->width;
	return R->laplacian[(scale * row) * R->width + (scale * column)];
}
float BoxIntegral(Matrix *img, int row, int col, int rows, int cols)
{
	float *data = img->data;
	int r1 = min(row, img->height) - 1;
	int c1 = min(col, img->width) - 1;
	int r2 = min(row + rows, img->height) - 1;
	int c2 = min(col + cols, img->width) - 1;
	float A = 0.0f, B = 0.0f, C = 0.0f, D = 0.0f;
	if (r1 >= 0 && c1 >= 0) A = data[r1 * img->width + c1];
	if (r1 >= 0 && c2 >= 0) B = data[r1 * img->width + c2];
	if (r2 >= 0 && c1 >= 0) C = data[r2 * img->width + c1];
	if (r2 >= 0 && c2 >= 0) D = data[r2 * img->width + c2];
	return max(0.f, A - B - C + D);
}
void buildResponseLayer(FastHessian *fh, ResponseLayer *rl)
{
	float *responses;
	unsigned char *laplacian;
	int step, b, l, w;
	float inverse_area, Dxx, Dyy, Dxy;
	int r, c, ar, index = 0, ac;
	responses = rl->responses;
	laplacian = rl->laplacian;
	step = rl->step;
	b = (rl->filter - 1) / 2;
	l = rl->filter / 3;
	w = rl->filter;
	inverse_area = 1.f / (float)(w*w);
	for (ar = 0; ar < rl->height; ++ar)
	{
		for (ac = 0; ac < rl->width; ++ac, index++)
		{
			r = ar * step;
			c = ac * step;
			Dxx = BoxIntegral(fh->img, r - l + 1, c - b, 2 * l - 1, w) - BoxIntegral(fh->img, r - l + 1, c - l / 2, 2 * l - 1, l) * 3.0f;
			Dyy = BoxIntegral(fh->img, r - b, c - l + 1, w, 2 * l - 1) - BoxIntegral(fh->img, r - l / 2, c - l + 1, l, 2 * l - 1) * 3.0f;
			Dxy = BoxIntegral(fh->img, r - l, c + 1, l, l) + BoxIntegral(fh->img, r + 1, c - l, l, l) - BoxIntegral(fh->img, r - l, c - l, l, l) - BoxIntegral(fh->img, r + 1, c + 1, l, l);
			Dxx *= inverse_area;
			Dyy *= inverse_area;
			Dxy *= inverse_area;
			responses[index] = (Dxx * Dyy - 0.81f * Dxy * Dxy);
			laplacian[index] = (Dxx + Dyy >= 0 ? 1 : 0);
		}
	}
}
void allocateResponseMap(FastHessian *fh)
{
	int w = (fh->i_width / fh->init_sample);
	int h = (fh->i_height / fh->init_sample);
	int s = fh->init_sample;
	if (fh->octaves >= 1)
	{
		ResLayer(&fh->responseMap[0], w, h, s, 9);
		ResLayer(&fh->responseMap[1], w, h, s, 15);
		ResLayer(&fh->responseMap[2], w, h, s, 21);
		ResLayer(&fh->responseMap[3], w, h, s, 27);
	}
	if (fh->octaves >= 2)
	{
		ResLayer(&fh->responseMap[4], w >> 1, h >> 1, s << 1, 39);
		ResLayer(&fh->responseMap[5], w >> 1, h >> 1, s << 1, 51);
	}
	if (fh->octaves >= 3)
	{
		ResLayer(&fh->responseMap[6], w >> 2, h >> 2, s << 2, 75);
		ResLayer(&fh->responseMap[7], w >> 2, h >> 2, s << 2, 99);
	}
	if (fh->octaves >= 4)
	{
		ResLayer(&fh->responseMap[8], w >> 3, h >> 3, s << 3, 147);
		ResLayer(&fh->responseMap[9], w >> 3, h >> 3, s << 3, 195);
	}
	if (fh->octaves >= 5)
	{
		ResLayer(&fh->responseMap[10], w >> 4, h >> 4, s << 4, 291);
		ResLayer(&fh->responseMap[11], w >> 4, h >> 4, s << 4, 387);
	}
}
void FastHessianInit(FastHessian *fhts, Matrix *img, const int octaves, const int init_sample, const float thres)
{
	fhts->img = img;
	fhts->i_height = img->height;
	fhts->i_width = img->width;
	fhts->octaves = octaves;
	fhts->init_sample = init_sample;
	fhts->thresh = thres;
	fhts->iptssize = 0;
	allocateResponseMap(fhts);
}
int isExtremum(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, FastHessian* fh)
{
	int rr, cc;
	int layerBorder = (t->filter + 1) / (t->step << 1);
	if (r <= layerBorder || r >= t->height - layerBorder || c <= layerBorder || c >= t->width - layerBorder)
		return 0;
	float candidate = getResponse(r, c, t, m);
	if (candidate < fh->thresh)
		return 0;
	for (rr = -1; rr <= 1; ++rr)
	{
		for (cc = -1; cc <= 1; ++cc)
		{
			if (get_Response(r + rr, c + cc, t) >= candidate || ((rr != 0 || cc != 0) && getResponse(r + rr, c + cc, t, m) >= candidate) || getResponse(r + rr, c + cc, t, b) >= candidate)
				return 0;
		}
	}
	return 1;
}
void deriv3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, double *dI)
{
	dI[0] = (getResponse(r, c + 1, t, m) - getResponse(r, c - 1, t, m)) / 2.0;
	dI[1] = (getResponse(r + 1, c, t, m) - getResponse(r - 1, c, t, m)) / 2.0;
	dI[2] = (get_Response(r, c, t) - getResponse(r, c, t, b)) / 2.0;
}
void hessian3D(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, double H[3][3])
{
	double v = getResponse(r, c, t, m);
	double dxx = getResponse(r, c + 1, t, m) + getResponse(r, c - 1, t, m) - 2.0 * v;
	double dyy = getResponse(r + 1, c, t, m) + getResponse(r - 1, c, t, m) - 2.0 * v;
	double dss = get_Response(r, c, t) + getResponse(r, c, t, b) - 2.0 * v;
	double dxy = (getResponse(r + 1, c + 1, t, m) - getResponse(r + 1, c - 1, t, m) - getResponse(r - 1, c + 1, t, m) + getResponse(r - 1, c - 1, t, m)) / 4.0;
	double dxs = (get_Response(r, c + 1, t) - get_Response(r, c - 1, t) - getResponse(r, c + 1, t, b) + getResponse(r, c - 1, t, b)) / 4.0;
	double dys = (get_Response(r + 1, c, t) - get_Response(r - 1, c, t) - getResponse(r + 1, c, t, b) + getResponse(r - 1, c, t, b)) / 4.0;
	H[0][0] = dxx;
	H[0][1] = dxy;
	H[0][2] = dxs;
	H[1][0] = dxy;
	H[1][1] = dyy;
	H[1][2] = dys;
	H[2][0] = dxs;
	H[2][1] = dys;
	H[2][2] = dss;
}
void solve3x3Fast(const double A[3][3], const double b[3], double x[3])
{
	double invdet = 1.0 / (DBL_EPSILON + (A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1])
		- A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0])
		+ A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0])));
	x[0] = invdet * (b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) -
		A[0][1] * (b[1] * A[2][2] - A[1][2] * b[2]) +
		A[0][2] * (b[1] * A[2][1] - A[1][1] * b[2]));
	x[1] = invdet * (A[0][0] * (b[1] * A[2][2] - A[1][2] * b[2]) -
		b[0] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) +
		A[0][2] * (A[1][0] * b[2] - b[1] * A[2][0]));
	x[2] = invdet * (A[0][0] * (A[1][1] * b[2] - b[1] * A[2][1]) -
		A[0][1] * (A[1][0] * b[2] - b[1] * A[2][0]) +
		b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]));
}
void interpolateExtremumFast(int r, int c, ResponseLayer *t, ResponseLayer *m, ResponseLayer *b, FastHessian *fh)
{
	double dI[3], H[3][3], out[3];
	deriv3D(r, c, t, m, b, dI);
	hessian3D(r, c, t, m, b, H);
	solve3x3Fast(H, dI, out);
	out[0] = -out[0];
	out[1] = -out[1];
	out[2] = -out[2];
	if (fabs(out[0]) < 0.5 && fabs(out[1]) < 0.5 && fabs(out[2]) < 0.5)
	{
		if (!fh->iptssize)
		{
			fh->iptssize = 384;
			fh->ipts = (Ipoint*)malloc(sizeof(Ipoint) * fh->iptssize);
		}
		if (fh->iptssize == fh->interestPtsLen)
		{
			fh->iptssize <<= 1;
			fh->ipts = (Ipoint*)realloc(fh->ipts, sizeof(Ipoint) * fh->iptssize);
		}
		fh->ipts[fh->interestPtsLen].x = ((float)c + (float)out[0]) * (float)t->step;
		fh->ipts[fh->interestPtsLen].y = ((float)r + (float)out[1]) * (float)t->step;
		fh->ipts[fh->interestPtsLen].scale = 0.1333f * ((float)m->filter + (float)out[2] * (float)(m->filter - b->filter));
		fh->ipts[fh->interestPtsLen].laplacian = getLaplacian(r, c, t, m);
		fh->interestPtsLen++;
	}
}
void getIpoints(FastHessian *fh)
{
	ResponseLayer *b, *m, *t;
	static const int filter_map[5][4] = { { 0,1,2,3 },{ 1,3,4,5 },{ 3,5,6,7 },{ 5,7,8,9 },{ 7,9,10,11 } };
	for (int o = 0; o < fh->octaves; ++o)
		for (int i = 0; i <= 1; ++i)
		{
			b = &fh->responseMap[filter_map[o][i]];
			m = &fh->responseMap[filter_map[o][i + 1]];
			t = &fh->responseMap[filter_map[o][i + 2]];
			for (int r = 0; r < t->height; ++r)
				for (int c = 0; c < t->width; ++c)
					if (isExtremum(r, c, t, m, b, fh))
						interpolateExtremumFast(r, c, t, m, b, fh);
		}
}
const float gauss25[7][7] = {// lookup table for 2d gaussian (sigma = 2.5) where (0,0) is top left and (6,6) is bottom right
0.02546481f, 0.02350698f, 0.01849125f, 0.01239505f, 0.00708017f, 0.00344629f, 0.00142946f,
0.02350698f, 0.02169968f, 0.01706957f, 0.01144208f, 0.00653582f, 0.00318132f, 0.00131956f,
0.01849125f, 0.01706957f, 0.01342740f, 0.00900066f, 0.00514126f, 0.00250252f, 0.00103800f,
0.01239505f, 0.01144208f, 0.00900066f, 0.00603332f, 0.00344629f, 0.00167749f, 0.00069579f,
0.00708017f, 0.00653582f, 0.00514126f, 0.00344629f, 0.00196855f, 0.00095820f, 0.00039744f,
0.00344629f, 0.00318132f, 0.00250252f, 0.00167749f, 0.00095820f, 0.00046640f, 0.00019346f,
0.00142946f, 0.00131956f, 0.00103800f, 0.00069579f, 0.00039744f, 0.00019346f, 0.00008024f };
float haarX(Matrix *img, const int row, const int column, const int s, const int sD)
{
	return BoxIntegral(img, row - sD, column, s, sD) - 1.0f * BoxIntegral(img, row - sD, column - sD, s, sD);
}
float haarY(Matrix *img, const int row, const int column, const int s, const int sD)
{
	return BoxIntegral(img, row, column - sD, sD, s) - 1.0f * BoxIntegral(img, row - sD, column - sD, sD, s);
}
#ifndef M_PI
#define M_PI 3.141592653589793f
#endif
#ifndef M_2PI
#define M_2PI 6.283185307179864f
#endif
float getAngle(float X, float Y)
{
	if (X > 0.0f && Y >= 0.0f)
		return atanf(Y / X);
	if (X < 0.0f && Y >= 0.0f)
		return M_PI - atanf(-Y / X);
	if (X < 0.0f && Y < 0.0f)
		return M_PI + atanf(Y / X);
	if (X > 0.0f && Y < 0.0f)
		return M_2PI - atanf(-Y / X);
	return 0;
}
#define ANGRESLEN 109
void getOrientation(SURFDescriptor *despt, int index)
{
	float gauss = 0.f, scale = despt->hessians->ipts[index].scale;
	const int s = fRound(scale), r = fRound(despt->hessians->ipts[index].y), c = fRound(despt->hessians->ipts[index].x);
	float resX[ANGRESLEN] = { 0.0f }, resY[ANGRESLEN] = { 0.0f }, Ang[ANGRESLEN] = { 0.0f };
	const int id[] = { 6,5,4,3,2,1,0,1,2,3,4,5,6 };
	int idx = 0;
	for (int i = -6; i <= 6; ++i)
	{
		for (int j = -6; j <= 6; ++j)
		{
			if (i*i + j * j < 36)
			{
				gauss = gauss25[id[i + 6]][id[j + 6]];
				resX[idx] = gauss * haarX(despt->integralImg, r + j * s, c + i * s, s << 2, s << 1);
				resY[idx] = gauss * haarY(despt->integralImg, r + j * s, c + i * s, s << 2, s << 1);
				Ang[idx] = getAngle(resX[idx], resY[idx]);
				++idx;
			}
		}
	}
	float sumX, sumY;
	float max = 0.f, orientation = 0.f;
	float ang1, ang2 = 0.f;
	for (ang1 = 0.0f; ang1 < M_2PI; ang1 += 0.15f)
	{
		ang2 = (ang1 + 1.04719755f > M_2PI ? ang1 - 5.23598775f : ang1 + 1.04719755f); //ang2 = (ang1 + M_PI / 3.0f > M_2PI ? ang1 - 5.0f*M_PI / 3.0f : ang1 + M_PI / 3.0f);
		sumX = sumY = 0.f;
		for (unsigned int k = 0; k < ANGRESLEN; ++k)
		{
			const float ang = Ang[k];
			if (ang1 < ang2 && ang1 < ang && ang < ang2)
			{
				sumX += resX[k];
				sumY += resY[k];
			}
			else if (ang2 < ang1 && ((ang > 0 && ang < ang2) || (ang > ang1 && ang < M_2PI)))
			{
				sumX += resX[k];
				sumY += resY[k];
			}
		}
		if (sumX*sumX + sumY * sumY > max)
		{
			max = sumX * sumX + sumY * sumY;
			orientation = getAngle(sumX, sumY);
		}
	}
	despt->hessians->ipts[index].orientation = orientation;
}
float gaussian1(int x, int y, float sig)
{
	return (1.0f / (M_2PI*sig*sig)) * expf(-(x*x + y * y) / (2.0f*sig*sig));
}
float gaussian2(float x, float y, float sig)
{
	return 1.0f / (M_2PI*sig*sig) * expf(-(x*x + y * y) / (2.0f*sig*sig));
}
void getDescriptor(SURFDescriptor *despt, int index)
{
	int y, x, sample_x, sample_y, count = 0, roundedScale;
	int i, ix = 0, j = 0, jx = 0, xs = 0, ys = 0;
	float scale, *desc, dx, dy, mdx, mdy, co, si;
	float gauss_s1 = 0.f, gauss_s2 = 0.f;
	float rx = 0.f, ry = 0.f, rrx = 0.f, rry = 0.f, len = 0.f;
	float cx = -0.5f, cy = 0.f;
	scale = despt->hessians->ipts[index].scale;
	roundedScale = fRound(scale);
	x = fRound(despt->hessians->ipts[index].x);
	y = fRound(despt->hessians->ipts[index].y);
	desc = despt->hessians->ipts[index].descriptor;
	co = cosf(despt->hessians->ipts[index].orientation);
	si = sinf(despt->hessians->ipts[index].orientation);
	i = -8;
	while (i < 12)
	{
		j = -8;
		i = i - 4;
		cx += 1.f;
		cy = -0.5f;
		while (j < 12)
		{
			dx = dy = mdx = mdy = 0.f;
			cy += 1.f;
			j = j - 4;
			ix = i + 5;
			jx = j + 5;
			xs = fRound(x + (-jx * scale*si + ix * scale*co));
			ys = fRound(y + (jx*scale*co + ix * scale*si));
			for (int k = i; k < i + 9; ++k)
			{
				for (int l = j; l < j + 9; ++l)
				{
					sample_x = fRound(x + (-l * scale*si + k * scale*co));
					sample_y = fRound(y + (l*scale*co + k * scale*si));
					gauss_s1 = gaussian1(xs - sample_x, ys - sample_y, 2.5f * scale);
					rx = haarX(despt->integralImg, sample_y, sample_x, roundedScale << 1, roundedScale);
					ry = haarY(despt->integralImg, sample_y, sample_x, roundedScale << 1, roundedScale);
					rrx = gauss_s1 * (-rx * si + ry * co);
					rry = gauss_s1 * (rx*co + ry * si);
					dx += rrx;
					dy += rry;
					mdx += fabsf(rrx);
					mdy += fabsf(rry);
				}
			}
			gauss_s2 = gaussian2(cx - 2.0f, cy - 2.0f, 1.5f);
			desc[count++] = dx * gauss_s2;
			desc[count++] = dy * gauss_s2;
			desc[count++] = mdx * gauss_s2;
			desc[count++] = mdy * gauss_s2;
			len += (dx*dx + dy * dy + mdx * mdx + mdy * mdy) * gauss_s2*gauss_s2;
			j += 9;
		}
		i += 9;
	}
	len = 1.0f / sqrtf(len);
	for (int i = 0; i < 64; ++i)
		desc[i] *= len;
}
typedef struct
{
	int threadPerTask, offset;
	SURFDescriptor *surf;
} genSURFThread;
void *threadResLayer(void *args)
{
	genSURFThread *arguments = (genSURFThread*)args;
	for (int i = 0; i < arguments->threadPerTask; ++i)
	{
		getOrientation(arguments->surf, i + arguments->offset);
		getDescriptor(arguments->surf, i + arguments->offset);
	}
	return 0;
}
void getDescriptors(SURFDescriptor *despt)
{
	if (!despt->hessians->interestPtsLen)
		return;
	int spawnThreads = 3;
	if (spawnThreads > 3)
		spawnThreads = 3;
	pthread_t th[3];
	genSURFThread thData[3];
	int threadPerTask = despt->hessians->interestPtsLen / (spawnThreads + 1);
	thData[0].offset = threadPerTask + (despt->hessians->interestPtsLen - threadPerTask * (spawnThreads + 1));
	thData[0].surf = despt;
	thData[0].threadPerTask = threadPerTask;
	for (int j = 1; j < spawnThreads; j++)
	{
		thData[j].offset = thData[0].offset + threadPerTask * j;
		thData[j].surf = despt;
		thData[j].threadPerTask = threadPerTask;
	}
	for (int j = 0; j < spawnThreads; j++)
		pthread_create(&th[j], 0, threadResLayer, (void*)&thData[j]);
	for (int i = 0; i < thData[0].offset; ++i)
	{
		getOrientation(despt, i);
		getDescriptor(despt, i);
	}
	for (int j = 0; j < spawnThreads; j++)
		pthread_join(th[j], 0);
}
void IpointPair_add(IpointPair *v, Ipoint *first, Ipoint *second)
{
	size_t memSize;
	if (!v->size)
	{
		v->size = 10;
		memSize = sizeof(Ipoint) * v->size;
		v->pt1 = (Ipoint*)malloc(memSize);
		v->pt2 = (Ipoint*)malloc(memSize);
		memset(v->pt1, 0, memSize);
		memset(v->pt2, 0, memSize);
	}
	if (v->size == v->count)
	{
		v->size <<= 1;
		memSize = sizeof(Ipoint) * v->size;
		v->pt1 = (Ipoint*)realloc(v->pt1, memSize);
		v->pt2 = (Ipoint*)realloc(v->pt2, memSize);
	}
	memcpy(&v->pt1[v->count], first, sizeof(Ipoint));
	memcpy(&v->pt2[v->count], second, sizeof(Ipoint));
	v->count++;
}
IpointPair getMatches(SURFDescriptor *surf1, SURFDescriptor *surf2, float threshold)
{
	IpointPair pair = { 0 };
	float dist, d1, d2;
	Ipoint *match;
	int i, j, idx;
	for (i = 0; i < surf1->hessians->interestPtsLen; i++)
	{
		d1 = d2 = FLT_MAX;
		for (j = 0; j < surf2->hessians->interestPtsLen; j++)
		{
			dist = 0.0f;
			for (idx = 0; idx < 64; ++idx)
				dist += (surf1->hessians->ipts[i].descriptor[idx] - surf2->hessians->ipts[j].descriptor[idx])*(surf1->hessians->ipts[i].descriptor[idx] - surf2->hessians->ipts[j].descriptor[idx]);
			dist = sqrtf(dist);
			if (dist < d1)
			{
				d2 = d1;
				d1 = dist;
				match = &surf2->hessians->ipts[j];
			}
			else if (dist < d2)
				d2 = dist;
		}
		if (d1 / d2 < threshold)
		{
			surf1->hessians->ipts[i].dx = match->x - surf1->hessians->ipts[i].x;
			surf1->hessians->ipts[i].dy = match->y - surf1->hessians->ipts[i].y;
			IpointPair_add(&pair, &surf1->hessians->ipts[i], match);
		}
	}
	return pair;
}
void matcherFree(IpointPair *pair)
{
	free(pair->pt1);
	free(pair->pt2);
}
void Integral(float *img, Matrix *integralImg)
{
	float rs = 0.0f;
	for (int j = 0; j < integralImg->width; j++)
	{
		rs += img[j];
		integralImg->data[j] = rs;
	}
	for (int i = 1; i < integralImg->height; ++i)
	{
		rs = 0.0f;
		for (int j = 0; j < integralImg->width; ++j)
		{
			rs += img[i*integralImg->width + j];
			integralImg->data[i*integralImg->width + j] = rs + integralImg->data[(i - 1)*integralImg->width + j];
		}
	}
}
void SURFInitialize(SURFDescriptor *surf, Matrix *img, int octaves, int init_sample, float thres)
{
	surf->hessians = (FastHessian*)malloc(sizeof(FastHessian));
	surf->integralImg = MatrixAllocatePtr(img->width, img->height);
	memset(surf->hessians->responseMap, 0, 12 * sizeof(ResponseLayer));
	FastHessianInit(surf->hessians, surf->integralImg, octaves, init_sample, thres);
	surf->hessians->ipts = 0;
}
void SURFCompute(SURFDescriptor *surf, float *img)
{
	Integral(img, surf->integralImg);
	for (int i = 0; i < 12; ++i)
		buildResponseLayer(surf->hessians, &surf->hessians->responseMap[i]);
	if (surf->hessians->interestPtsLen > 0)
		free(surf->hessians->ipts);
	surf->hessians->iptssize = 0;
	surf->hessians->interestPtsLen = 0;
	// Extract interest points and store in vector ipts
	getIpoints(surf->hessians);
	getDescriptors(surf);
}
void SURFFree(SURFDescriptor *surf)
{
	free(surf->integralImg->data);
	free(surf->integralImg);
	if (surf->hessians->interestPtsLen > 0)
		free(surf->hessians->ipts);
	for (int i = 0; i < 12; i++)
	{
		if (surf->hessians->responseMap[i].responses)
		{
			free(surf->hessians->responseMap[i].responses);
			free(surf->hessians->responseMap[i].laplacian);
		}
	}
	free(surf->hessians);
}