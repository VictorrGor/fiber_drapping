

//--------------------------------------------------------------------------------------
struct VS_INPUT
{
	float4 Pos		: SV_POSITION;
	float4 WorldPos : POSITION;
	float4 Color	: COLOR;
	float3 Normal	: NORMAL;
};



struct PointLight 
{ 
	float4 Ambient; 
	float4 Diffuse; 
	float4 Specular; 
	float3 Position; 
	float Range; 
	float3 Att; 
	float pad; 
};

struct Material
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular;
	float4 Reflect;
};

struct light
{
	float3 dir;
	float pad;
	float4 ambient;
	float4 diffuse;
};

cbuffer cbPerFrame : register(b0)
{
	//PointLight gPointLight;
	light ll;
	//float3 gEyePosW;
	//float4x4 mVWP;
};

cbuffer cbPerObject : register(b1)
{
	float4 Ambient;
	float4 Diffuse;
	float4 Specular;
	float4 Reflect;
	bool isMaterial[16];
}

//https://habr.com/ru/post/441862/
//https://ru.wikipedia.org/wiki/%D0%97%D0%B0%D1%82%D0%B5%D0%BD%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BF%D0%BE_%D0%A4%D0%BE%D0%BD%D0%B3%D1%83
float4 main(VS_INPUT input) : SV_Target
{
	if (isMaterial[0])
	{
		input.Normal = normalize(input.Normal);
		float4 diffuse = Diffuse * ll.ambient;
		float3 res = diffuse +saturate(dot(ll.dir, input.Normal) * ll.diffuse * Diffuse);
		return float4(res, Diffuse.a);
	//	float4 plPos = mul(gPointLight.Position, mVWP);

	//	float3 finalColor = float3(0.0f, 0.0f, 0.0f);

	//	float3 lightVec = plPos - input.Pos.xyz;
	//	float3 toEye = gEyePosW - input.Pos.xyz;

	//	float4 finalAmbient = mat.Ambient * gPointLight.Ambient;
	//	float len = length(lightVec);
	//	
	//	len = length(input.Pos.xyz);
	//	if (len > 10) return float4(1,0,0,1);
	//	if (input.Pos.y > 1) return float4(0, 1, 0, 1);
	//	if (input.Pos.z > 1) return float4(0, 0, 1, 1);
	//	if (len > 100) return float4(0.9, 0.1, 0.1, 1);//finalAmbient;
	//	lightVec /= len;

	//	
	//	float4 res = mat.Ambient * float4(0, 0, 0, 0) ; /*Фонового освещениея нет*/
	//	float ll = dot(lightVec, input.Normal);
	//	res = res + mat.Diffuse * ll;//diffuse
	//	res += mat.Specular * dot(-lightVec, toEye);
	//	return float4(finalColor, mat.Diffuse.a);
	}
	return input.Color;
}
