#pragma once

struct Color{
    unsigned char b, g, r;
};

Color operator+(const Color& c1, const Color& c2)
{
    auto adjust_and_check = [](const unsigned char& value1, const unsigned char& value2){
        int temp = value1 + value2;
        if(temp < 0)
            temp = 0;
        else if(temp > 255)
            temp = 255;
        
        return static_cast<unsigned char>(temp);
    };

    return Color{adjust_and_check(c1.b, c2.b), adjust_and_check(c1.g, c2.g), adjust_and_check(c1.r, c2.r)};
}

Color& operator*(double d, Color& c)
{
    c.b *= d;
    c.g *= d;
    c.r *= d;

    return c;
}