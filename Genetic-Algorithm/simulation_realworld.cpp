#include <algorithm>
#include <cmath>
#include <iostream>

#include "simulation.h"

int util::h(int t)
{
    if (t >= constant::T_w_0 && t <= constant::T_w_0 + 24)
    {
        return std::min(constant::T_w, t - constant::T_w_0);
    }
    else if (t < constant::T_w_0)
    {
        return h(t + 24) - constant::T_w;
    }
    else
        return h(t - 24) + constant::T_w;
}

int util::g(int t1, int t2)
{
    return util::h(t1 + t2) - util::h(t2);
}

util::matrix<double, K, K> util::D2T(util::matrix<int, K, K> D, int speed)
{
    util::matrix<double, K, K> T;

    for (size_t i = 0; i < K; i++)
        for (size_t j = 0; j < K; j++)
            T[i][j] = std::ceil(D[i][j] / static_cast<double>(speed));

    return T;
}

void sim::simulation(const int data[M][K], double (&time)[S], const int T_0, const int ship_speed)
{

    /*
        이 코드는 임시적으로 거리 배열과 시나리오가 달라진 코드를 반영하기 위한 코드입니다.
    */
    using util::D2T;
    using util::g;
    using util::h;

    int tick_count;

    const util::matrix<int, K, K> distance = {
        {
            {19, 25, 13, 123, 173, 229, 301, 242, 290, 322, 362, 359, 392, 441, 522, 569},
            {31, 39, 14, 112, 162, 219, 291, 232, 280, 311, 351, 348, 381, 431, 511, 558},
            {55, 54, 29, 88, 138, 195, 267, 208, 256, 288, 328, 325, 357, 407, 488, 535},
            {81, 80, 55, 62, 112, 169, 241, 182, 230, 262, 301, 299, 331, 381, 461, 509},
            {103, 102, 77, 32, 93, 150, 221, 163, 211, 242, 282, 279, 312, 362, 442, 489},
            {120, 119, 95, 56, 73, 130, 201, 143, 191, 222, 262, 259, 292, 342, 422, 469},
            {176, 175, 151, 112, 46, 108, 182, 70, 70, 207, 233, 230, 262, 312, 393, 440},
            {235, 233, 209, 170, 99, 30, 14, 38, 86, 129, 154, 152, 184, 234, 315, 362},

        }};

    auto T = D2T(distance, ship_speed);

    const std::array<int, S> NEW_SCENARIO = {29600, 29600, 29600, 27500, 27500, 20700, 27500, 21700};

    // auto T = D2T(constant::D, ship_speed);

    // 각 지역 사이의 이동거리와 배의 속도를 고려한 이동 시간
    // 각 지역에서 사건 발생 지점까지 이동하는데 걸리는 시간

    /*
    std::cout << "distance: " << std::endl;
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
            std::cout << distance[i][j] << ' ';
        std::cout << std::endl;
    }
    std::cout << "T: " << std::endl;
    for (int i = 0; i < K; i++)
    {
        for (int j = 0; j < K; j++)
            std::cout << T[i][j] << ' ';
        std::cout << std::endl;
    }
    std::cout << "NEW_SCENARIO: " << std::endl;
    for (int i = 0; i < S; i++)
        std::cout << NEW_SCENARIO[i] << std::endl;
    */
    for (size_t i = 0; i < S; i++) // S: 16개의 지역을 대상으로 진행
    {
        tick_count = 0; // tick_count는 일한 시간을 의미한다.
        double Q_i = NEW_SCENARIO[i];
        // double Q_i = SCENARIO[i];

        if (i >= 8)
            break;

        auto target = Q_i / (constant::alpha * constant::beta * constant::gamma) / 3.0;
        std::array<int, K> work_times;
        for (size_t j = 0; j < K; j++)
        {
            work_times[j] = -g(T[i][j], T_0 + T[i][i]);
        }

        int day1_work_acc = 0;
        int day2_work_acc = 0;

        int first = true;
        int second = true;
        while (true)
        {
            tick_count++;
            int work_acc = 0;
            /*
            std::cout << std::endl;
            std::cout << "1. tick_count : " << tick_count << std::endl;
            */
            for (size_t j = 0; j < K; j++)
            {
                work_times[j]++;

                if (tick_count <= 8) // 첫 날을 의미한다.
                {
                    // work_acc += (std::max(work_times[j], 0) * data[0][j] * 1); // 첫 날의 회수 효율은 12/7
                    work_acc += (std::max(work_times[j], 0) * data[0][j] * 12 / 7);
                    if (tick_count == 8 && j == K - 1)
                    {
                        day1_work_acc = work_acc;
                        // 8시간 동안 즉, 하루 동안 마지막 지역의 노동력까지 다 이용했을 때, 회수량을 기록한다.
                        for (size_t k = 0; k < K; k++)
                        {
                            if (work_times[k] > 1)
                                work_times[k] = 0;
                            /*
                            예를 들어, 8시간(하루)동안 각 지역의 인력들이 노동한 시간이 4, 5, 8, 0, 0, 0, -5, 0, -8, -8, -4, -8, -8, -8, -13, -16 이라면
                            노동한 시간, 회수 효율, 각 지역의 장비의 수를 곱해서 회수량을 저장한다.다음 날은 회수 효율이 떨어지기 때문에
                            이 회수량을 day1_work_acc에 기록한 후 위 노동 시간을 0, 0, 0, 0, 0, 0, -5, 0, -8, -8, -4, -8, -8, -8, -13, -16으로 변환
                            그럼 다음 날은 1, 1, 1, 1, 1, 1, -4, 1, -7, -7, -3, -7, -7, -7, -12, -15로 변환 된 후 회수량을 계산하고 첫 날 회수량을 더하게 된다.
                            */
                        }
                    }
                }

                else if (tick_count <= 16) // 둘 째날을 의미한다.
                {
                    // work_acc += (std::max(work_times[j], 0) * data[0][j] * 1);
                    work_acc += (std::max(work_times[j], 0) * data[0][j] * 6 / 7);
                    if (first)
                    {
                        first = false;
                        work_acc += day1_work_acc;
                    }
                    if (j == K - 1)
                        first = true;
                    if (tick_count == 16 && j == K - 1)
                    {
                        day2_work_acc = work_acc;
                        // 16시간 동안 즉, 이틀 동안 마지막 지역의 노동력까지 다 이용했을 때, 회수량을 기록한다.
                        for (size_t k = 0; k < K; k++)
                        {
                            if (work_times[k] > 1)
                                work_times[k] = 0;
                        }
                    }
                }
                else // 셋 째날을 의미한다.
                {
                    // work_acc += (std::max(work_times[j], 0) * data[0][j] * 1);
                    work_acc += (std::max(work_times[j], 0) * data[0][j] * 3 / 7);
                    if (second)
                    {
                        second = false;
                        work_acc += day2_work_acc;
                    }
                    if (j == K - 1)
                        second = true;
                }
            }

            /*
            std::cout << "2. work_rate " << work_acc << "/" << target << std::endl;
            std::cout << "3. work_times: ";
            for (size_t j = 0; j < K; j++)
            {
                std::cout << work_times[j] << " ";
            }
            std::cout << std::endl;
            */
            if (target <= work_acc) // 회수량이 목표량보다 크거나 같을 경우, 함수 종료
            {
                time[i] = std::ceil(tick_count * constant::gamma) / 8.0;
                break;
            }
        }
    }

    return;
}

int main()
{
    /* testcase1
    int data[1][K] = {1268, 1082, 1363, 1360, 867, 315, 2622, 514, 723, 917, 788, 2134, 2958, 736, 748, 245};
    int data[1][K] = {1772, 250, 9375, 792, 1772, 125, 9375, 167, 105, 355, 250, 522, 9375, 167, 105, 11};
    int data[1][K] = {3403, 50, 3230, 1457, 1967, 2424, 4431, 179, 50, 284, 50, 371, 4282, 222, 50, 50};
    int data[1][K] = {1714, 109, 8012, 22, 1131, 11, 6726, 14, 10, 245, 121, 363, 4000, 11, 10, 1};
    int data[1][K] = {1320, 114, 8477, 249, 861, 11, 6842, 11, 10, 309, 238, 11, 4000, 36, 10, 1};
    int data[1][K] = {1422, 11, 8325, 11, 1166, 11, 6897, 11, 10, 277, 11, 251, 4075, 11, 10, 1};
    int data[1][K] = {1616, 244, 754, 769, 396, 121, 7727, 153, 83, 335, 237, 492, 9358, 162, 51, 2};
    */
    /* testcase2 (수정 전)
    int data[1][K] = {2833, 0, 3400, 1267, 2833, 200, 4200, 267, 0, 533, 0, 833, 5700, 267, 167, 0}; // 기존 배치

    int data[1][K] = {713, 184, 5534, 133, 548, 11, 4000, 11, 10, 105, 187, 45, 4000, 11, 10, 1}; // GA 배치 (1:1:1) (8, 10)
    int data[1][K] = {1340, 264, 8964, 589, 11, 174, 4281, 237, 138, 297, 369, 600, 4375, 245, 10, 1}; // GA 배치 (4:2:1) (8, 10)
    int data[1][K] = {2691, 181, 9775, 1129, 2757, 149, 5682, 37, 24, 252, 187, 517, 4812, 213, 20, 16}; // GA 배치 (1:1:1) (8, 5)
    int data[1][K] = {1630, 36, 9819, 716, 1459, 96, 8087, 54, 66, 63, 222, 374, 5619, 62, 10, 1}; // GA 배치 (4:2:1) (8, 5)
    int data[1][K] = {11, 11, 4959, 11, 417, 11, 4000, 11, 10, 11, 11, 11, 4000, 11, 10, 1}; // GA 배치 (1:1:1) (12, 10)
    int data[1][K] = {1066, 11, 5271, 11, 704, 11, 4000, 11, 10, 11, 11, 11, 4000, 11, 10, 1}; // GA 배치 (4:2:1) (12, 10)
    int data[1][K] = {1911, 344, 9677, 1209, 2682, 152, 4706, 209, 124, 188, 335, 556, 4915, 225, 120, 1}; // GA 배치 (1:1:1) (12, 5)
    int data[1][K] = {1003, 174, 9924, 784, 1454, 105, 7508, 147, 85, 235, 230, 508, 5343, 111, 96, 2}; // GA 배치 (4:2:1) (12, 5)
    */

    /*
    int data[1][K] = {755, 109, 5899, 188, 11, 24, 4000, 11, 10, 11, 129, 115, 4158, 11, 10, 1}; // GA 배치 (1:1:1) (8, 10)
    int data[1][K] = {2617, 328, 9855, 996, 2580, 138, 5886, 62, 28, 11, 58, 292, 5394, 11, 40, 1}; // GA 배치 (1:1:1) (8, 5)
    int data[1][K] = {11, 221, 5123, 11, 11, 11, 4000, 11, 10, 11, 11, 11, 4000, 11, 10, 1}; // GA 배치 (1:1:1) (12, 10)
    int data[1][K] = {2172, 280, 9835, 976, 2756, 200, 4092, 266, 112, 316, 400, 653, 4905, 266, 78, 16}; // GA 배치 (1:1:1) (12, 5)
    int data[1][K] = {1255, 237, 6000, 778, 1246, 125, 4920, 59, 50, 350, 231, 508, 4827, 130, 23, 11}; // GA 배치 (4:2:1) (8, 10)
    int data[1][K] = {1503, 119, 9963, 789, 1514, 125, 7514, 61, 105, 338, 250, 372, 6074, 119, 18, 11}; // GA 배치 (4:2:1) (8, 5)
    int data[1][K] = {272, 29, 6040, 775, 856, 11, 4000, 20, 10, 11, 11, 37, 4094, 11, 10, 1}; // GA 배치 (4:2:1) (12, 10)
    int data[1][K] = {1720, 171, 8915, 735, 1711, 75, 7947, 118, 43, 286, 208, 467, 5571, 121, 54, 1}; // GA 배치 (4:2:1) (12, 5)
    */
    double time[S];
    int data[1][K] = {272, 29, 6040, 775, 856, 11, 4000, 20, 10, 11, 11, 37, 4094, 11, 10, 1}; // GA 배치 (4:2:1) (12, 10)
    sim::simulation(data, time, 12, 10);
    /*
        시작 시간과 직원들이 일 할 수 있는 시간, 방제 자원을 사용할 수 있는 시간을 고려하여 소요 시간을 구한다.
        작업일수는 3일, 작업시간은 하루 당 8시간이다.
        방제자원 동원 시간은 해당지역 장비의 출항 준비 시간과 항해시간
        인근지역 장비의 출항 준비 시간과 항해시간

        지역의 거리를 고려하여, 코드를 수정
        int 형 배열 dist를 입력으로 받는다.
   */
    for (int i = 0; i < K; ++i)
    {
        std::cout << time[i] * 8 << ',';
    }
}
