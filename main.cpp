
#define WIN32_LEAN_AND_MEAN

#define MIN(x, y) ((x) <= (y) ? (x) : (y))
#define MAX(x, y) ((x) >= (y) ? (x) : (y))

#include <Windows.h>

#include "stdint.h"

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <filesystem>
#include <functional>
#include <list>
#include <map>
#include <memory>
#include <mutex>
#include <queue>
#include <random>
#include <regex>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <combaseapi.h>

#include <dxgi1_6.h>
#include <d3d11_4.h>
#include <d3d11shader.h>
#include <d3dcompiler.h>
#include <directxmath.h>

#pragma comment(lib, "d3d11.lib")
#pragma comment(lib, "dwmapi.lib")
#pragma comment(lib, "dxguid.lib")

using namespace std;

namespace FastLED
{
	typedef uint8_t   fract8;

	static uint8_t const p[] = {
	151, 160, 137,  91,  90,  15, 131,  13, 201,  95,  96,  53, 194, 233,   7, 225,
	140,  36, 103,  30,  69, 142,   8,  99,  37, 240,  21,  10,  23, 190,   6, 148,
	247, 120, 234,  75,   0,  26, 197,  62,  94, 252, 219, 203, 117,  35,  11,  32,
	 57, 177,  33,  88, 237, 149,  56,  87, 174,  20, 125, 136, 171, 168,  68, 175,
	 74, 165,  71, 134, 139,  48,  27, 166,  77, 146, 158, 231,  83, 111, 229, 122,
	 60, 211, 133, 230, 220, 105,  92,  41,  55,  46, 245,  40, 244, 102, 143,  54,
	 65,  25,  63, 161,   1, 216,  80,  73, 209,  76, 132, 187, 208,  89,  18, 169,
	200, 196, 135, 130, 116, 188, 159,  86, 164, 100, 109, 198, 173, 186,   3,  64,
	 52, 217, 226, 250, 124, 123,   5, 202,  38, 147, 118, 126, 255,  82,  85, 212,
	207, 206,  59, 227,  47,  16,  58,  17, 182, 189,  28,  42, 223, 183, 170, 213,
	119, 248, 152,   2,  44, 154, 163,  70, 221, 153, 101, 155, 167,  43, 172,   9,
	129,  22,  39, 253,  19,  98, 108, 110,  79, 113, 224, 232, 178, 185, 112, 104,
	218, 246,  97, 228, 251,  34, 242, 193, 238, 210, 144,  12, 191, 179, 162, 241,
	 81,  51, 145, 235, 249,  14, 239, 107,  49, 192, 214,  31, 181, 199, 106, 157,
	184,  84, 204, 176, 115, 121,  50,  45, 127,   4, 150, 254, 138, 236, 205,  93,
	222, 114,  67,  29,  24,  72, 243, 141, 128, 195,  78,  66, 215,  61, 156, 180,
	151 };

#define P(x) (*((const  uint8_t*)(p + (x))))

	int8_t avg7(int8_t i, int8_t j)
	{
		return (i >> 1) + (j >> 1) + (i & 0x1);
	}

	uint8_t scale8(uint8_t i, fract8 scale)
	{
#if 1 // (FASTLED_SCALE8_FIXED == 1)
		return (((uint16_t)i) * (1 + (uint16_t)(scale))) >> 8;
#else
		return ((uint16_t)i * (uint16_t)(scale)) >> 8;
#endif
	}

	uint8_t ease8InOutQuad(uint8_t i)
	{
		uint8_t j = i;
		if (j & 0x80) {
			j = 255 - j;
		}
		uint8_t jj = scale8(j, j);
		uint8_t jj2 = jj << 1;
		if (i & 0x80) {
			jj2 = 255 - jj2;
		}
		return jj2;
	}
#define EASE8(x)  (ease8InOutQuad(x) )

	int8_t grad8(uint8_t hash, int8_t x)
	{
		// since the tests below can be done bit-wise on the bottom
		// four bits, there's no need to mask off the higher bits
		//	hash = hash & 15;

		int8_t u, v;
		if (hash & 8) {
			u = x; v = x;
		}
		else {
			if (hash & 4) {
				u = 1; v = x;
			}
			else {
				u = x; v = 1;
			}
		}

		if (hash & 1) { u = -u; }
		if (hash & 2) { v = -v; }

		return avg7(u, v);
	}
	int8_t  grad8(uint8_t hash, int8_t x, int8_t y)
	{
		// since the tests below can be done bit-wise on the bottom
		// three bits, there's no need to mask off the higher bits
		//  hash = hash & 7;

		int8_t u, v;
		if (hash & 4) {
			u = y; v = x;
		}
		else {
			u = x; v = y;
		}

		if (hash & 1) { u = -u; }
		if (hash & 2) { v = -v; }

		return avg7(u, v);
	}

	uint8_t qadd8(uint8_t i, uint8_t j)
	{
		unsigned int t = i + j;
		if (t > 255) t = 255;
		return t;
	}


	int8_t inline  lerp7by8(int8_t a, int8_t b, fract8 frac)
	{
		// int8_t delta = b - a;
		// int16_t prod = (uint16_t)delta * (uint16_t)frac;
		// int8_t scaled = prod >> 8;
		// int8_t result = a + scaled;
		// return result;
		int8_t result;
		if (b > a) {
			uint8_t delta = b - a;
			uint8_t scaled = scale8(delta, frac);
			result = a + scaled;
		}
		else {
			uint8_t delta = a - b;
			uint8_t scaled = scale8(delta, frac);
			result = a - scaled;
		}
		return result;
	}
	int8_t inoise8_raw(uint16_t x)
	{
		// Find the unit cube containing the point
		uint8_t X = x >> 8;

		// Hash cube corner coordinates
		uint8_t A = P(X);
		uint8_t AA = P(A);
		uint8_t B = P(X + 1);
		uint8_t BA = P(B);

		// Get the relative position of the point in the cube
		uint8_t u = (uint8_t)x;

		// Get a signed version of the above for the grad function
		int8_t xx = ((uint8_t)(x) >> 1) & 0x7F;
		uint8_t N = 0x80;

		u = EASE8(u);

		int8_t ans = lerp7by8(grad8(P(AA), xx), grad8(P(BA), xx - N), u);

		return ans;
	}

	int8_t inoise8_raw(uint16_t x, uint16_t y)
	{
		// Find the unit cube containing the point
		uint8_t X = x >> 8;
		uint8_t Y = y >> 8;

		// Hash cube corner coordinates
		uint8_t A = P(X) + Y;
		uint8_t AA = P(A);
		uint8_t AB = P(A + 1);
		uint8_t B = P(X + 1) + Y;
		uint8_t BA = P(B);
		uint8_t BB = P(B + 1);

		// Get the relative position of the point in the cube
		uint8_t u = (uint8_t)x;
		uint8_t v = (uint8_t)y;

		// Get a signed version of the above for the grad function
		int8_t xx = ((uint8_t)(x) >> 1) & 0x7F;
		int8_t yy = ((uint8_t)(y) >> 1) & 0x7F;
		uint8_t N = 0x80;

		u = EASE8(u); v = EASE8(v);

		int8_t X1 = lerp7by8(grad8(P(AA), xx, yy), grad8(P(BA), xx - N, yy), u);
		int8_t X2 = lerp7by8(grad8(P(AB), xx, yy - N), grad8(P(BB), xx - N, yy - N), u);

		int8_t ans = lerp7by8(X1, X2, v);

		return ans;
		// return scale8((70+(ans)),234)<<1;
	}

	uint8_t inoise8(uint16_t x, uint16_t y) {
		//return scale8(69+inoise8_raw(x,y),237)<<1;
		int8_t n = inoise8_raw(x, y);  // -64..+64
		n += 64;                         //   0..128
		uint8_t ans = qadd8(n, n);     //   0..255
		return ans;
	}

	void fill_raw_noise8(uint8_t *pData, uint8_t num_points, uint8_t octaves, uint16_t x, int scale, uint16_t time) {
		uint32_t _xx = x;
		uint32_t scx = scale;
		for (int o = 0; o < octaves; ++o) {
			for (int i = 0, xx = _xx; i < num_points; ++i, xx += scx) {
				pData[i] = qadd8(pData[i], inoise8(xx, time) >> o);
			}

			_xx <<= 1;
			scx <<= 1;
		}
	}
}

namespace
{

	inline uint8_t min3(uint8_t a, uint8_t b, uint8_t c)
	{
		return MIN(a, MIN(b, c));
	}

	inline void rgbw_2_rgb(const uint8_t *rgbw, uint8_t *rgb)
	{
		rgb[0] = MAX(rgbw[0], rgbw[3]);
		rgb[1] = MAX(rgbw[1], rgbw[3]);
		rgb[2] = MAX(rgbw[2], rgbw[3]);
	}

	inline void rgb_2_rgbw(const uint8_t *rgb, uint8_t *rgbw)
	{
		uint8_t min_component = min3(rgb[0], rgb[1], rgb[2]);
		if (min_component <= 84)
		{
			rgbw[0] = rgb[0] - min_component;
			rgbw[1] = rgb[1] - min_component;
			rgbw[2] = rgb[2] - min_component;
			rgbw[3] = 3 * min_component;
		}
		else
		{
			uint8_t w3 = 85;
			rgbw[0] = rgb[0] - w3;
			rgbw[1] = rgb[1] - w3;
			rgbw[2] = rgb[2] - w3;
			rgbw[3] = 255;
		}
	}

	class Timer
	{
	public:
		Timer()
		{
			QueryPerformanceFrequency(&freq);
			QueryPerformanceCounter(&t);
		}

		inline float Update()
		{
			LARGE_INTEGER t2;
			QueryPerformanceCounter(&t2);
			float dt = (t2.QuadPart - t.QuadPart) / float(freq.QuadPart);
			t = t2;
			return dt;
		}

		//inline float UpdateFromReferenceTime(LARGE_INTEGER &tExt)
		//{
		//	float dt = (tExt.QuadPart - t.QuadPart) / float(freq.QuadPart);
		//	t = tExt;
		//	return dt;
		//}

		inline void UpdateReferenceTime()
		{
			QueryPerformanceCounter(&t);
		}

		LARGE_INTEGER t;
		LARGE_INTEGER freq;
	};


	void copySmall(uint8_t *dst, uint8_t *src, int size)
	{
		--size;
		for (; size >= 0; --size)
		{
			dst[size] = src[size];
		}
	}

	uint64_t xorShift64(uint64_t a)
	{
		a ^= (a << 21);
		a ^= (a >> 35);
		a ^= (a << 4);
		return a;
	}

	uint64_t state = xorShift64(0xCAFEBABE);

	uint64_t nextLong()
	{
		uint64_t a = state;
		state = xorShift64(a);
		return a;
	}


	D3D11_VIEWPORT vp = {};
	ID3D11Device *device;
	ID3D11DeviceContext *context;
	IDXGISwapChain *swapChain;
	ID3D11RenderTargetView *rtv = nullptr;
	ID3D11ShaderResourceView *rsv = nullptr;
	ID3D11Texture2D *ledTexture;
	D3D11_TEXTURE2D_DESC texDescLast = {};
	uint8_t *rgba = nullptr;
	Timer timer;
	float totalTime = 0.0;

#define TEX_PADDING 4

#define LED_BUFFER_SIZE_SIDE 109
#define LED_BUFFER_SIZE_BACK 251
#define LED_BUFFER_SIZE_TOTAL (LED_BUFFER_SIZE_SIDE * 2 + LED_BUFFER_SIZE_BACK)
	uint8_t ledBuffer[LED_BUFFER_SIZE_TOTAL * 4] = {};
	uint8_t *ledBufferCollar = ledBuffer;

	LRESULT CALLBACK WindowProcedure(HWND hWnd, UINT message, WPARAM wParam, LPARAM lParam)
	{
		switch (message)
		{
		case WM_TIMER:
		{
			if (wParam == 100)
			{
				//RECT rect;
				//GetClientRect(hWnd, &rect);
				//vp.Width = (float)rect.right - rect.left;
				//vp.Height = (float)rect.bottom - rect.top;

				float dt = timer.Update();

				context->RSSetViewports(1, &vp);

				float ClearColor[4] = { 0.0f, 0.0f, 0.0f, 1.0f };
				context->ClearRenderTargetView(rtv, ClearColor);
				context->OMSetRenderTargets(1, &rtv, nullptr);


				D3D11_TEXTURE2D_DESC finalTexDesc;
				ledTexture->GetDesc(&finalTexDesc);

				if (finalTexDesc.Width != texDescLast.Width ||
					finalTexDesc.Height != texDescLast.Height)
				{
					delete[] rgba;
					rgba = new uint8_t[finalTexDesc.Width * finalTexDesc.Height * 4];
					texDescLast = finalTexDesc;
				}

				memset(rgba, 0, finalTexDesc.Width * finalTexDesc.Height * 4);

				// Do animation

				totalTime += dt;

				uint8_t noiseNumbers[LED_BUFFER_SIZE_TOTAL] = {};
				uint16_t offset = uint16_t(totalTime * 500.0f);
				for (int i = 0; i < LED_BUFFER_SIZE_TOTAL; i += 255)
				{
					int numLeft = MIN(255, LED_BUFFER_SIZE_TOTAL - i);
					FastLED::fill_raw_noise8(noiseNumbers + i, numLeft, 1, offset + i * 96, 96, uint16_t(totalTime * 0.0f));
				}

				//for (int i = 0; i < LED_BUFFER_SIZE_TOTAL; ++i)
				//{
				//	noiseNumbers[i] = noiseNumbers[i] >> 1;
				//}

				for (int i = 0; i < LED_BUFFER_SIZE_TOTAL; ++i)
				{
					ledBufferCollar[i * 4 + 0] = noiseNumbers[i];
					ledBufferCollar[i * 4 + 1] = noiseNumbers[i];
					ledBufferCollar[i * 4 + 2] = noiseNumbers[i];
					ledBufferCollar[i * 4 + 3] = noiseNumbers[i];
				}

				// Apply led buffer to texture buffer

				// Right side
				for (int i = 0; i < LED_BUFFER_SIZE_SIDE; ++i)
				{
					uint8_t *read = ledBuffer + i * 4;
					int wx = TEX_PADDING + LED_BUFFER_SIZE_BACK - 1;
					int wy = TEX_PADDING + LED_BUFFER_SIZE_SIDE - i;
					uint8_t *write = rgba + (wy * texDescLast.Width + wx) * 4;
					rgbw_2_rgb(read, write);
				}

				// Left side
				for (int i = 0; i < LED_BUFFER_SIZE_SIDE; ++i)
				{
					uint8_t *read = ledBuffer + (LED_BUFFER_SIZE_SIDE + LED_BUFFER_SIZE_BACK + LED_BUFFER_SIZE_SIDE - i - 1) * 4;
					int wx = TEX_PADDING;
					int wy = TEX_PADDING + LED_BUFFER_SIZE_SIDE - i;
					uint8_t *write = rgba + (wy * texDescLast.Width + wx) * 4;
					rgbw_2_rgb(read, write);
				}

				// Back
				for (int i = 0; i < LED_BUFFER_SIZE_BACK; ++i)
				{
					uint8_t *read = ledBuffer + (LED_BUFFER_SIZE_SIDE + LED_BUFFER_SIZE_BACK - i - 1) * 4;
					int wx = TEX_PADDING + i;
					int wy = TEX_PADDING;
					uint8_t *write = rgba + (wy * texDescLast.Width + wx) * 4;
					rgbw_2_rgb(read, write);
				}

				context->UpdateSubresource(ledTexture, 0, 0, rgba, texDescLast.Width * 4, 0);
				context->Draw(6, 0);
				swapChain->Present(0, 0);
			}
		}
		return 0;

		case WM_DESTROY:
		{
			delete[] rgba;
			PostQuitMessage(0);
		}
		break;
		}

		return DefWindowProc(hWnd, message, wParam, lParam);
	}
}

int WINAPI WinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR szCmdLine, int iCmdShow)
{
	HRESULT hr;

	CoInitializeEx(nullptr, COINIT_APARTMENTTHREADED | COINIT_DISABLE_OLE1DDE);

	WNDCLASSEX windowClass = { 0 };
	windowClass.cbSize = sizeof(WNDCLASSEX);
	windowClass.hInstance = GetModuleHandle(nullptr);
	windowClass.hbrBackground = (HBRUSH)(COLOR_BACKGROUND);
	windowClass.lpfnWndProc = WindowProcedure;
	windowClass.lpszClassName = L"TestClass";
	windowClass.style = CS_PARENTDC;
	auto wndClass = RegisterClassEx(&windowClass);

	RECT rect = { 200, 200, 3000, 1600 };
	DWORD exStyle = 0;
	HWND hWnd = CreateWindowEx(exStyle, (LPCTSTR)wndClass, L"LEDTest",
		WS_OVERLAPPEDWINDOW | WS_VISIBLE,
		rect.left, rect.top, rect.right - rect.left, rect.bottom - rect.top,
		nullptr, 0, GetModuleHandle(nullptr), 0);

	GetClientRect(hWnd, &rect);
	vp.Width = (float)rect.right - rect.left;
	vp.Height = (float)rect.bottom - rect.top;

	UINT createDeviceFlags = D3D11_CREATE_DEVICE_DEBUG;

	const D3D_FEATURE_LEVEL featureLevels[] =
	{
		D3D_FEATURE_LEVEL_11_1,
	};
	const UINT numFeatureLevels = ARRAYSIZE(featureLevels);


	D3D_FEATURE_LEVEL acceptedLevel;
	hr = D3D11CreateDevice(nullptr, D3D_DRIVER_TYPE_HARDWARE, nullptr, createDeviceFlags, featureLevels, numFeatureLevels,
		D3D11_SDK_VERSION, &device, &acceptedLevel, &context);

	IDXGIFactory1 *dxgiFactory = nullptr;
	IDXGIDevice *dxgiDevice = nullptr;

	hr = device->QueryInterface(__uuidof(IDXGIDevice), reinterpret_cast<void**>(&dxgiDevice));
	if (SUCCEEDED(hr))
	{
		IDXGIAdapter *adapter = nullptr;
		hr = dxgiDevice->GetAdapter(&adapter);
		if (SUCCEEDED(hr))
		{
			hr = adapter->GetParent(__uuidof(IDXGIFactory1), reinterpret_cast<void**>(&dxgiFactory));
			adapter->Release();
		}
		dxgiDevice->Release();
	}

	const bool useFlipPresentation = true;

	DXGI_SWAP_CHAIN_DESC sd = { 0 };
	sd.BufferCount = useFlipPresentation ? 2 : 1;
	sd.BufferDesc.Width = 0;
	sd.BufferDesc.Height = 0;
	sd.BufferDesc.Format = DXGI_FORMAT_R8G8B8A8_UNORM;
	sd.BufferUsage = DXGI_USAGE_RENDER_TARGET_OUTPUT;
	sd.OutputWindow = hWnd;
	sd.Windowed = TRUE;
	sd.SwapEffect = useFlipPresentation ? DXGI_SWAP_EFFECT_FLIP_DISCARD : DXGI_SWAP_EFFECT_DISCARD;
	sd.SampleDesc.Count = 1;
	sd.SampleDesc.Quality = 0;
	hr = dxgiFactory->CreateSwapChain(device, &sd, &swapChain);
	dxgiFactory->Release();

	ID3D11Texture2D *backBufferTexture;
	if (SUCCEEDED(swapChain->GetBuffer(0, __uuidof(ID3D11Texture2D), reinterpret_cast<void**>(&backBufferTexture))))
	{
		device->CreateRenderTargetView(backBufferTexture, nullptr, &rtv);
		backBufferTexture->Release();
	}

	D3D11_TEXTURE2D_DESC texDesc = { 0 };
	texDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
	texDesc.SampleDesc.Count = 1;
	texDesc.ArraySize = 1;
	texDesc.MipLevels = 1;
	texDesc.Usage = D3D11_USAGE_DEFAULT;
	texDesc.CPUAccessFlags = D3D11_CPU_ACCESS_WRITE;
	texDesc.BindFlags = D3D11_BIND_SHADER_RESOURCE;
	texDesc.Width = TEX_PADDING * 2 + LED_BUFFER_SIZE_BACK;
	texDesc.Height = TEX_PADDING * 2 + 1 + LED_BUFFER_SIZE_SIDE;
	hr = device->CreateTexture2D(&texDesc, nullptr, &ledTexture);

	D3D11_SHADER_RESOURCE_VIEW_DESC resDesc;
	resDesc.Format = DXGI_FORMAT_B8G8R8A8_UNORM;
	resDesc.ViewDimension = D3D_SRV_DIMENSION_TEXTURE2D;
	resDesc.Texture2D.MipLevels = 1;
	resDesc.Texture2D.MostDetailedMip = 0;
	hr = device->CreateShaderResourceView(ledTexture, &resDesc, &rsv);





	ID3D11SamplerState *samplerState;
	ID3D11Buffer *vertexBuffer;
	ID3D11VertexShader *vertexShader;
	ID3D11PixelShader *pixelShader;
	ID3D11InputLayout *inputLayout;


	D3D11_SAMPLER_DESC samplerDesc;

	const int vertexCount = 6;
	UINT strideInBytes = (3 + 2) * 4;
	UINT offset = 0;
	const float vData[vertexCount * (3 + 2)] = {
		-1, 1, 0, 0, 0,
		1, -1, 0, 1, 1,
		-1, -1, 0, 0, 1,

		1, 1, 0, 1, 0,
		1, -1, 0, 1, 1,
		-1, 1, 0, 0, 0,
	};

	D3D11_BUFFER_DESC bd = { 0 };
	D3D11_SUBRESOURCE_DATA initData = { 0 };
	D3D11_BLEND_DESC blendState = { 0 };
	D3D11_DEPTH_STENCIL_DESC depthStencilDesc = { 0 };
	D3D11_RASTERIZER_DESC rasterizerStateDesc = { D3D11_FILL_SOLID, D3D11_CULL_NONE };

	memset(&samplerDesc, 0, sizeof(samplerDesc));
	samplerDesc.Filter = D3D11_FILTER_MIN_MAG_MIP_POINT; //D3D11_FILTER_MIN_MAG_MIP_LINEAR;
	samplerDesc.AddressU = samplerDesc.AddressV = samplerDesc.AddressW = D3D11_TEXTURE_ADDRESS_CLAMP;
	samplerDesc.MaxAnisotropy = 1;
	samplerDesc.ComparisonFunc = D3D11_COMPARISON_ALWAYS;
	samplerDesc.MaxLOD = D3D11_FLOAT32_MAX;
	hr = device->CreateSamplerState(&samplerDesc, &samplerState);


	bd.Usage = D3D11_USAGE_IMMUTABLE;
	bd.ByteWidth = strideInBytes * vertexCount;
	bd.BindFlags = D3D11_BIND_VERTEX_BUFFER;
	initData.pSysMem = vData;
	hr = device->CreateBuffer(&bd, &initData, &vertexBuffer);




	HMODULE moduleCompiler = LoadLibrary(L"d3dcompiler_47.dll");

	pD3DCompile fnCompile = (pD3DCompile)GetProcAddress(moduleCompiler, "D3DCompile");

	ID3DBlob *blobVS = nullptr;
	ID3DBlob *blobPS = nullptr;
	ID3DBlob *errorMessages = nullptr;

	vector<D3D11_INPUT_ELEMENT_DESC> layout =
	{
		{ "POSITION", 0, DXGI_FORMAT_R32G32B32_FLOAT, 0,  0, D3D11_INPUT_PER_VERTEX_DATA, 0 },
		{ "TEXCOORD", 0, DXGI_FORMAT_R32G32_FLOAT,    0, 12, D3D11_INPUT_PER_VERTEX_DATA, 0 },
	};

	string srcVS = R"(struct VS_INPUT
{
	float3 a_Position : POSITION;
	float2 a_TexCoord : TEXCOORD0;
};

struct VS_OUTPUT
{
	float4 v_Position : SV_POSITION;
	float2 v_TexCoord : TEXCOORD0;
};

VS_OUTPUT main(VS_INPUT IN)
{
	VS_OUTPUT OUT;

	OUT.v_Position = float4(IN.a_Position, 1.0);
	OUT.v_TexCoord = IN.a_TexCoord;
	return OUT;
})";
	string srcPS = R"(struct VS_OUTPUT
{
	float4 v_Position : SV_POSITION;
	float2 v_TexCoord : TEXCOORD0;
};

struct PS_OUTPUT
{
	float4 p_FragColor : SV_TARGET;
};

Texture2D g_Texture0:register(t0);
SamplerState g_Texture0SamplerState:register(s0);

PS_OUTPUT main(VS_OUTPUT IN)
{
	PS_OUTPUT OUT;
	OUT.p_FragColor = g_Texture0.SampleLevel(g_Texture0SamplerState, IN.v_TexCoord, 0.0); //float4(0.3, 0.0, 0.3, 1.0);
	return OUT;
})";


	hr = fnCompile(srcVS.c_str(), srcVS.size(), "t_vs", nullptr, nullptr, "main", "vs_4_0",
		D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &blobVS, &errorMessages);
	if (errorMessages)
		errorMessages->Release();
	errorMessages = nullptr;

	hr = fnCompile(srcPS.c_str(), srcPS.size(), "t_ps", nullptr, nullptr, "main", "ps_4_0",
		D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3, 0, &blobPS, &errorMessages);
	if (errorMessages)
		errorMessages->Release();
	errorMessages = nullptr;

	hr = device->CreateVertexShader(blobVS->GetBufferPointer(),
		blobVS->GetBufferSize(), nullptr, &vertexShader);

	hr = device->CreatePixelShader(blobPS->GetBufferPointer(),
		blobPS->GetBufferSize(), nullptr, &pixelShader);

	hr = device->CreateInputLayout(layout.data(), (UINT)layout.size(), blobVS->GetBufferPointer(),
		blobVS->GetBufferSize(), &inputLayout);


	blobPS->Release();
	blobVS->Release();


	SetTimer(hWnd, 100, 33, nullptr);

	//shared->ctx->RSSetState(shared->rasterizerState);
	//shared->ctx->OMSetBlendState(shared->blendState, nullptr, 0xFFFFFFFF);

	context->IASetInputLayout(inputLayout);
	context->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);

	//context->VSSetConstantBuffers(0, 1, &shaderConstantBuffer);
	//context->PSSetConstantBuffers(0, 1, &shaderConstantBuffer);

	context->VSSetShader(vertexShader, nullptr, 0);
	context->PSSetShader(pixelShader, nullptr, 0);

	ID3D11SamplerState *samplers[1] = { samplerState };
	ID3D11ShaderResourceView *resourceViews[1] = { rsv };

	context->PSSetSamplers(0, 1, samplers);
	context->PSSetShaderResources(0, 1, resourceViews);

	context->IASetVertexBuffers(0, 1, &vertexBuffer, &strideInBytes, &offset);

	timer.UpdateReferenceTime();

	MSG msg = { 0 };
	while (GetMessage(&msg, 0, 0, 0))
	{
		TranslateMessage(&msg);
		DispatchMessage(&msg);
	}

	inputLayout->Release();
	pixelShader->Release();
	vertexShader->Release();
	samplerState->Release();
	vertexBuffer->Release();
	rsv->Release();
	ledTexture->Release();
	rtv->Release();
	swapChain->Release();
	context->Release();
	device->Release();

	FreeLibrary(moduleCompiler);

	CoUninitialize();

	return 0;
}
