// This is a personal academic project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com
#include "Input/Mouse.h"

Mouse::Mouse()
{
	isLeftKeyPressed = false;
	isRightKeyPressed = false;
	wheel_pos = 1;
	mousePos = POINT();
	savedMPos = POINT();
}

void Mouse::updWheelPos(int newPos)
{
	wheel_pos = newPos;
}

void Mouse::updLK(bool isPressed)
{
	isLeftKeyPressed = isPressed;
	savedMPos = mousePos;
}

void Mouse::updRK(bool isPressed)
{
	isRightKeyPressed = isPressed;
}

void Mouse::updMousePos(POINT mPos)
{
	if (!isLeftKeyPressed)
	{
		savedMPos = mPos;
	}
	mousePos = mPos;
}

POINT Mouse::updateSavedPos()
{
	POINT res;
	res.x = mousePos.x - savedMPos.x;
	res.y = mousePos.y - savedMPos.y;
	savedMPos = mousePos;
	return res;
}

int Mouse::getWheelPos()
{
	return wheel_pos;
}
