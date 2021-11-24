// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#define DEBUG_CONSOLE
#include "Render/Render.h"

#include "Render/Utils.h"

RenderSys rs;

void incrementRenderObjCount()
{
	rs.incrementRenderObjCount();
}
void decrementRenderObjCount()
{
	rs.decrementRenderObjCount();
}

LRESULT CALLBACK WndProc(HWND, UINT, WPARAM, LPARAM);

HRESULT InitWindow(HINSTANCE hInst, int nCmdShow, HWND* hWnd)
{
	convertTextInfoFileToBin("resources/texture/ASCI_desc.txt", "resources/texture/ASCI_desc.binary");
	file.open("log", std::ios_base::trunc);
	WNDCLASSEX wc = { 0 };
	wchar_t name[] = L"Name";

	wc.cbSize = sizeof(wc);
	wc.hInstance = hInst;
	wc.style = CS_HREDRAW | CS_VREDRAW;
	wc.lpfnWndProc = WndProc;
	wc.cbClsExtra = 0;
	wc.cbWndExtra = 0;
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);
	wc.hCursor = LoadCursor(NULL, IDC_ARROW);
	wc.hbrBackground = (HBRUSH)GetStockObject(WHITE_BRUSH);
	wc.lpszMenuName = NULL;
	wc.lpszClassName = name;

	if (!RegisterClassEx(&wc))
	{
		MessageBox(NULL, L"Cann't register class", L"Error", MB_OK);
		return E_FAIL;
	}

	(*hWnd) = CreateWindow(name, L"Hello World", WS_OVERLAPPEDWINDOW, CW_USEDEFAULT, 0, CW_USEDEFAULT, 0, (HWND)NULL, (HMENU)NULL, (HINSTANCE)hInst, NULL);

	if (!hWnd)
	{
		MessageBox(NULL, L"Cann't create window", L"Error", MB_OK);
		return E_FAIL;
	}
	return S_OK;
}


void createMousePosBanner(POINT mousePt)
{
	char* str = new char[100];
	static TexturedObject* ar[2];
	static bool initFlag = true;
	Interface* pInterface = rs.getInterface();
	
	if (initFlag)
	{
		initFlag = false;
	}
	else
	{
		pInterface->deleteObject(ar[0]);
		pInterface->deleteObject(ar[1]);
	}
	snprintf(str, 100, "X: %d", mousePt.x);
	ar[0] = pInterface->makeWord(str, -1, 1, 4, rs.getDevice());

	snprintf(str, 100, "Y: %d", mousePt.y);
	ar[1] = pInterface->makeWord(str, -1, 0.9, 4, rs.getDevice());
	pInterface->pushTexturedObject(ar[0]);
	pInterface->pushTexturedObject(ar[1]);
	delete[] str;
}

void createObjectRenderCounter()
{
	char* str = new char[100];
	Interface* pInterface = rs.getInterface();
	static bool initFlag = true;
	static Button plus, minus;

	static TexturedObject* ar;
	if (initFlag)
	{
		initFlag = false;
		float plusX, plusY, minusX, minusY;
		plusX = -1;
		plusY = 0.7;
		TexturedObject* buf = pInterface->makeWord("+", plusX, plusY, 8, rs.getDevice());
		plus.setObj(buf);
		plus.setPayload(incrementRenderObjCount);
		///@todo Bad usage: bound box set arbitrarily. Window width and Height set by public fields of RenderSys, which must be private
		plus.setCoordinates(plusX, plusY, 50, 50, rs.width, rs.height);
		pInterface->pushActiveElement(&plus);

		minusX = -0.8;
		minusY = 0.7;
		buf = pInterface->makeWord("-", minusX, minusY, 8, rs.getDevice());
		pInterface->pushActiveElement(&minus);
		minus.setObj(buf);
		minus.setPayload(decrementRenderObjCount);
		minus.setCoordinates(minusX, minusY, 50, 50, rs.width, rs.height);
	}
	else
	{
		pInterface->deleteObject(ar);
	}
	snprintf(str, 100, "Objects count: %d", rs.getRenderCount());
	ar = pInterface->makeWord(str, -1, 0.8, 4, rs.getDevice());
	pInterface->pushTexturedObject(ar);

	delete[] str;
}

int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpCmdLine, int nCmdShow)
{
#ifdef DEBUG_CONSOLE
	if (!AllocConsole()) return 0;
	FILE* console = freopen("CONOUT$", "w", stdout);
#endif

	HWND hWnd;
	MSG msg = {0};

	if (FAILED(InitWindow(hInst, nCmdShow, &hWnd)))
	{
		return 0;
	}

	RECT rc;
	if (!GetClientRect(hWnd, &rc))
	{
		MessageBox(NULL, L"Getting rect faild! Something goes wrong.", L"Error!", MB_OK);
		rs.CleanupDevice();
		return 0;
	}
	UINT width = rc.right - rc.left;
	UINT height = rc.bottom - rc.top;
	rs.OnResize(width, height);

	if (FAILED(rs.InitDevice(&hWnd)))
	{
		MessageBox(NULL, L"Initialization faild! Something goes wrong.", L"Error!", MB_OK);
		rs.CleanupDevice();
		return 0;
	}

	ShowWindow(hWnd, nCmdShow);

	try
	{
		while (WM_QUIT != msg.message)
		{
			if (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE))
			{
				TranslateMessage(&msg);
				DispatchMessage(&msg);
			}
			else
			{
				rs.Render();
			}
		}
	}
	catch (...)
	{
		file.close();
	}

	rs.CleanupDevice();

	fclose(console);
	return msg.wParam;
}


LRESULT CALLBACK WndProc(HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
	HDC hDC;
	PAINTSTRUCT ps;

	static int wheel_dist = 1;
	
	switch (uMsg)
	{
		case WM_PAINT:
		{
			createObjectRenderCounter();
			hDC = BeginPaint(hWnd, &ps);
			EndPaint(hWnd, &ps);
			break;
		}
		case WM_CLOSE:
		{
			DestroyWindow(hWnd);
			break;
		}
		case WM_DESTROY:
		{
			file.close();
			PostQuitMessage(0);
			break;
		}
		case WM_KEYDOWN:
		{
			switch (wParam)
			{
				case VK_ESCAPE:
				{
					DestroyWindow(hWnd);
					break;
				}
				case VK_SPACE:
				{
					rs.disableMaterials();
					break;
				}
			};
			break;
		}
		case WM_MOUSEWHEEL:
		{
			int buff = (short)HIWORD(wParam);
			wheel_dist += 3 * (short)HIWORD(wParam) / WHEEL_DELTA;
			if (wheel_dist < 1) wheel_dist = 1;
			rs.getMouse()->updWheelPos(wheel_dist);
			break;
		}
		case WM_LBUTTONDOWN:
		{
			rs.getMouse()->updLK(true);
			break;
		}
		case WM_RBUTTONDOWN:
		{
			rs.getMouse()->updRK(true);
			break;
		}
		case WM_LBUTTONUP:
		{
			POINT pt;
			pt.x = LOWORD(lParam);
			pt.y = HIWORD(lParam);
			rs.getInterface()->SendClick(pt);
			rs.getMouse()->updLK(false);
			createObjectRenderCounter();
			break;
		}
		case WM_RBUTTONUP:
		{
			rs.getMouse()->updRK(false);
			break;
		}
		case WM_MOUSEMOVE:
		{
			POINT pt;
			pt.x = LOWORD(lParam);
			pt.y = HIWORD(lParam);
			rs.getMouse()->updMousePos(pt); 
			createMousePosBanner(pt);
			break;
		}
		case WM_SIZE:
		{
			UINT width = LOWORD(lParam);
			UINT height = HIWORD(lParam);
			rs.OnResize(width, height);
		}
		default:
			return DefWindowProc(hWnd, uMsg, wParam, lParam);
		}
	return 0;
}
