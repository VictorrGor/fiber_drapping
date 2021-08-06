#pragma once
#include <Windows.h>
#include <stdio.h>
#include <iostream>
#include "DataStructures.h"

//Convert text file with texture description to binary
//Text file format: 
//UINT count of characters
//char\t Left_u_float\t Right_u_float\t UINT_width

void deleteDoubleSpaces(char* str)
{
	int x = 0;
	for (int i = 0; str[i + 1]; ++i)
	{
		if ((!isspace(str[i])) || ((x > 0) && (!isspace(str[i + 1]))))
		{
			str[x] = str[i];
			++x;
		}
	}
	str[x] = '\0';
}

void convertTextInfoFileToBin(const char* filename_txt, const char* filename_bin)
{
	FILE* binary = fopen(filename_bin, "wb");
	if (!binary)
	{
		std::cout << "[convertTextInfoFileToBin] Cann't create binary file!\n";
		return;
	}

	FILE* text = fopen(filename_txt, "r");
	if (!text)
	{
		std::cout << "[convertTextInfoFileToBin] Cann't open text file!\n";
		fclose(binary);
		return;
	}

	TextTextureInfo tti;
	char* textBuf = new char[512];
	char* searchBuf = new char[512];
	char* value, *valueCopy;
	UINT diff = 0;
	UINT charactersCount = 0;

	char delim = ' ';

	//@todo Delete multiple spaces
	while (!feof(text))
	{
		fgets(textBuf, 512, text);
		deleteDoubleSpaces(textBuf);

		value = strchr(textBuf, (int)delim); //Read char
		diff = value - textBuf;
		if (diff == 1)
		{
			strncpy(&tti.ch, textBuf, diff);
		}
		else
		{
			throw "Size of first char more then one!";
		}
		value += 1;
		valueCopy = value;

		value = strchr(value, (int)delim); //Read left u value
		diff = value - valueCopy;
		strncpy(searchBuf, valueCopy, diff);
		searchBuf[diff] = '\0';
		tti.leftU = atof(searchBuf);
		value += 1;
		valueCopy = value;

		value = strchr(value, (int)delim); //Read right u value
		diff = value - valueCopy;
		strncpy(searchBuf, valueCopy, diff);
		searchBuf[diff] = '\0';
		tti.rightU = atof(searchBuf);
		value += 1;
		valueCopy = value;

		value = strchr(valueCopy, (int)'\n'); //Read width
		if (!value)
		{
			value = strchr(valueCopy, (int)'\0');
			if (!value)
			{
				throw "Utils.h::convertTextInfoFileToBin: Can't find null terminated character!\n";
			}
		}
		diff = value - valueCopy;	
		strncpy(searchBuf, valueCopy, diff);
		searchBuf[diff] = '\0';
		tti.pixelWidth = atoi(searchBuf);
		valueCopy = value;

		fwrite(&tti, sizeof(TextTextureInfo), 1, binary);
		++charactersCount;
	}

	tti.ch = ' ';
	tti.leftU = 0;
	tti.rightU = 0;
	tti.pixelWidth = 6;
	fwrite(&tti, sizeof(tti), 1, binary);
	++charactersCount;
	fwrite(&charactersCount, sizeof(UINT), 1, binary);

	fclose(binary);
	fclose(text);
}