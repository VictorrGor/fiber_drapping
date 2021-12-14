#pragma once
#include <Windows.h>
#pragma pack(1)
//class RenderSys;

class Mouse
{
	bool isLeftKeyPressed;
	bool isRightKeyPressed;
	int wheel_pos;

	POINT mousePos;
	POINT savedMPos;
public:
	//friend RenderSys;

	Mouse();
	void updWheelPos(int newPos);
	void updLK(bool isPressed);
	void updRK(bool isPressed);
	void updMousePos(POINT mPos);
	POINT updateSavedPos();
	int getWheelPos();
};
