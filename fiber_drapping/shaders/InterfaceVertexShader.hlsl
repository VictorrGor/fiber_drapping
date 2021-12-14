struct VS_INPUT
{
	float3 Pos : POSITION;
	float4 Color : COLOR;
	float3 Normal : NORMAL;
	float2 Tex: TEXCOORD0;
};

struct VS_OUTPUT
{
	float4 Pos : SV_POSITION;
	float4 Color : COLOR;
	float3 Normal : NORMAL;
	float2 Tex: TEXCOORD;
};

VS_OUTPUT main(VS_INPUT inp)
{
	VS_OUTPUT res = (VS_OUTPUT)0;
	
	res.Pos    = float4(inp.Pos, 1);
	res.Color  = inp.Color;
	res.Normal = inp.Normal;
	res.Tex    = inp.Tex;
	return res;
}