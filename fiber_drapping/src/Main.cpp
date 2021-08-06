#define DEBUG_CONSOLE
#include "Render.h"

#include "Utils.h"

RenderSys rs;

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
	wc.hIcon = LoadIcon(NULL, IDI_APPLICATION);

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



int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrev, LPSTR lpCmdLine, int nCmdShow)
{
#ifdef DEBUG_CONSOLE
	if (!AllocConsole()) return 0;
	freopen("CONOUT$", "w", stdout);
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
			rs.getMouse()->updLK(false);
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
