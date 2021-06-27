

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
	PointLight gPointLight;
	float3 gEyePosW;
	float dummy;
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
		float3 lightVec = gPointLight.Position - input.WorldPos.xyz;
		float3 toEyeW = normalize(gEyePosW - input.WorldPos.xyz);
		float lightVecDistance = length(lightVec);
		lightVec /= lightVecDistance;

		input.Normal = normalize(input.Normal);
		float4 ambient = Ambient * gPointLight.Ambient;
		float4 spec = float4(0, 0, 0, 0);
		float4 diff = float4(0, 0, 0, 0);
		
		if (lightVecDistance <= gPointLight.Range)
		{
			diff = (max(dot(input.Normal, lightVec), 0)) * gPointLight.Diffuse * Diffuse;

			float3 v = reflect(-lightVec, input.Normal);
			float specFactor = pow(max(dot(v, toEyeW), 0.0f), 32);
			spec = specFactor * Specular * gPointLight.Specular;
		}

		float3 res = diff + ambient + spec;
		return  float4(res, 1.0f);
	}
	return input.Color;
}
