/* 
tmp / data[j].b = 
- a * Sin(dyt) * Cos(dxt) / M_PI * ((1 - ef) / dxyh - r * pow(M_E, -r * dxyh)) / dxyh
- a * Sin(dyt) * Sin(dxt) * Sin(dxt) * M_PI * (ef * r / dxyh - (1 - ef) / dxyh / dxyh + r * r * ef) / dxyh / 2
+ a * Sin(dyt) * Sin(dxt) * Sin(dxt) / M_PI / M_PI / M_PI * ((1 - ef) / dxyh - r * ef) / dxyh / dxyh / 2

by using the following subsitutions:

dxt = 2*M_PI*x

(1 - cos(dxt)) = dxk

(1 - cos(dyt)) = dyk

dxk * M_PI_SQ = dxh

dyk * M_PI_SQ = dyh

(dxh + dyh ) / 2 = dxyh

pow(M_E, -r * dxyh) = ef

sD->KASQR = r

sD->A = a*/



tmp -= data[j].b * (-sD->A * cos(2 * M_PI * dx) / M_PI * sin(2 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 +
                    (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2))) /
                    ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2) - sD->KASQR * pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2))) /
                            ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2)
                    - sD->A * sin(2 * M_PI * dx) / (M_PI * M_PI) * sin(2 * M_PI * dy) * (pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 +
                    (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2)) *
                        sD->KASQR * sin(2 * M_PI * dx) / M_PI * log(M_E) / ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 +
                            (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2) -
                            (0.1e1 - pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) /
                                2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2))) *
                        pow((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 +
                        (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2, -2) *
                        sin(2 * M_PI * dx) / M_PI + sD->KASQR * sD->KASQR * pow(M_E, -sD->KASQR *
                        ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) /
                            2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2)) *
                        sin(2 * M_PI * dx) / M_PI * log(M_E)) / ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) /
                            2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2) / 2 + sD->A *
                    pow(sin(2 * M_PI * dx), 2) * pow(M_PI, -0.3e1) * sin(2 * M_PI * dy) * ((0.1e1 - pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 +
                    (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2))) /
                        ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2) - sD->KASQR * pow(M_E, -sD->KASQR * ((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2))) *
                    pow((1 - cos(2 * M_PI * dx)) / (M_PI * M_PI) / 2 + (1 - cos(2 * M_PI * dy)) / (M_PI * M_PI) / 2, -2) / 2) * multiplier;