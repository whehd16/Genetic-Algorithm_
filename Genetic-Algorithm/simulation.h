//
//  simulation.h
//  simulation
//
//  Created by jazz4rabbit on 11/9/18.
//

#include "ga_parameters.h"
#include <array>

#ifndef SIMULATION_H
#define SIMULATION_H

namespace util
{
template <class T, size_t ROW, size_t COL>
using matrix = std::array<std::array<T, COL>, ROW>;

int h(int t);
int g(int t1, int t2);
matrix<double, K, K> D2T(const matrix<int, K, K> D, const int speed);

} // namespace util

namespace constant
{
const int T_0 = 12; // 사건 발생시간
const int T_w_0 = 8;
const int T_w = 8;

const double alpha = 0.2;
const double beta = 1 / 3.0;
const double gamma = 1;

/*
	const util::matrix<int, K, K> D = {
		{
			{0, 0, 0, 0, 24, 24, 24, 24, 24, 48, 48, 48, 48, 48, 48, 48 },
			{ 0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,48 ,48 ,48 ,48 ,48 ,48 ,48 },
			{ 0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,48 ,48 ,48 ,48 ,48 ,48 ,48 },
			{ 0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,48 ,48 ,48 ,48 ,48 ,48 ,48 },
			{ 24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,24 ,24 },
			{ 24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,24 ,24 },
			{ 24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,24 ,24 },
			{ 24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,24 ,24 },
			{ 24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,24 ,24 ,24 ,24 ,24 ,24 ,24 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
			{ 48 ,48 ,48 ,48 ,24 ,24 ,24 ,24 ,24 ,0  ,0  ,0  ,0  ,0  ,0  ,0 },
		}
	};
	*/

const util::matrix<int, K, K> D = {

    // 인천, 평택, 대산, 군산, 목포, 완도, 여수, 제주, 서귀포, 통영, 창원, 부산, 울산, 포항, 동해, 속초 */
    /* 이전 거리 매트릭스
    {
			{0  ,41 ,32 ,122,187,270,309,264,294,347,403,406,435,486,527,604},
			{41 ,0  ,23 ,119,184,267,306,261,335,344,400,403,432,483,550,601},
			{32 ,23 ,0  ,99 ,160,93 ,281,236,317,319,275,378,407,458,525,576},
			{122,119,99 ,0  ,107,190,245,184,226,270,326,329,358,407,471,524},
			{187,184,160,107,0  ,83 ,151,96 ,109,189,225,228,257,309,382,433},
			{270,267,93 ,190,83 ,0  ,79 ,51 ,85 ,99 ,139,136,265,317,392,448},
			{309,306,281,245,151,79 ,0  ,108,122,42 ,101,97 ,126,164,241,293},
			{264,261,236,184,96 ,51 ,108,0  ,57 ,127,167,170,196,244,321,372},
			{294,335,317,226,109,85 ,122,57 ,0  ,164,177,219,248,286,363,415},
			{347,344,319,270,189,99 ,42 ,127,164,0  ,27 ,63 ,92 ,142,220,264},
			{403,400,275,326,225,139,101,167,177,27 ,0  ,41 ,69 ,125,196,253},
			{406,403,378,329,228,136,97 ,170,219,63 ,41 ,0  ,37 ,89 ,164,220},
			{435,432,407,358,257,265,126,196,248,92 ,69 ,37 ,0  ,56 ,137,185},
			{486,483,458,407,309,317,164,244,286,142,125,89 ,56 ,0  ,92 ,148},
			{527,550,525,471,382,392,241,321,363,220,196,164,137,92 ,0  ,57 },
			{604,601,576,524,433,448,293,372,415,264,253,220,185,148,57 ,0  },
		}
    */
    {
        {0, 41, 32, 122, 187, 270, 309, 264, 294, 347, 403, 406, 435, 486, 527, 604},
        {41, 0, 23, 119, 184, 267, 306, 261, 345, 344, 400, 403, 432, 483, 550, 601},
        {32, 23, 0, 99, 160, 243, 281, 236, 325, 319, 375, 378, 407, 458, 525, 576},
        {122, 119, 99, 0, 107, 190, 245, 184, 226, 270, 326, 329, 358, 407, 471, 524},
        {187, 184, 160, 107, 0, 83, 151, 96, 109, 189, 225, 228, 257, 309, 382, 433},
        {270, 267, 243, 190, 83, 0, 79, 51, 108, 99, 139, 138, 175, 227, 302, 358},
        {309, 306, 281, 245, 151, 79, 0, 108, 122, 42, 101, 97, 126, 164, 241, 293},
        {264, 261, 236, 184, 96, 51, 108, 0, 57, 127, 167, 170, 196, 244, 321, 372},
        {294, 345, 325, 226, 109, 108, 122, 57, 0, 184, 177, 180, 205, 259, 335, 400},
        {347, 344, 319, 270, 189, 99, 42, 127, 184, 0, 27, 63, 92, 142, 220, 264},
        {403, 400, 375, 326, 225, 139, 101, 167, 177, 27, 0, 41, 69, 125, 195, 253},
        {406, 403, 378, 329, 228, 138, 97, 170, 180, 63, 41, 0, 37, 89, 164, 220},
        {435, 432, 407, 358, 257, 175, 126, 196, 205, 92, 69, 37, 0, 56, 137, 185},
        {486, 483, 458, 407, 309, 227, 164, 244, 259, 142, 125, 89, 56, 0, 92, 148},
        {527, 550, 525, 471, 382, 302, 241, 321, 335, 220, 195, 164, 137, 92, 0, 57},
        {604, 601, 576, 524, 433, 358, 293, 372, 400, 264, 253, 220, 185, 148, 57, 0},
    }}; // namespace constant

const int ship_speed = 10;

} // namespace constant

namespace sim
{
/*
	 *  0:인천,    1:평택,  2:대산,  3:군산,
	 *  4:목포,    5:완도,  6:여수,  7:제주,
	 *  8:서귀포,  9:통영, 10:창원, 11:부산,
	 * 12:울산,   13:포항, 14:동해, 15:속초

	constexpr int S = 16; // # of simulation
	const std::array<int, S> SCENARIO = {
		8500,  1200, 45000, 3800,
		8500,   600, 45000,  800,
		500,   1700,  1200, 2500,
		45000,  800,   500,   50
	};
*/
constexpr int S = 16; // # of simulation

const std::array<int, S> SCENARIO = {
    8500,
    1200,
    45000,
    3800,
    8500,
    600,
    45000,
    800,
    500,
    1700,
    1200,
    2500,
    45000,
    800,
    500,
    50

};

/*
	 * Maritime Equipment(해상 장비)
	 * Ground Equipment(육상 장비?)
	 * Oil-skimmer(유회수기)
	*/
constexpr int cur_mar_skimmer[S][1] = {
    569, 814, 594, 600,
    739, 220, 1667, 305,
    309, 287, 474, 1662,
    1617, 450, 505, 200};
constexpr int cur_gro_skimmer[S][1] = {
    699, 268, 769, 760,
    128, 95, 955, 209,
    414, 630, 314, 472,
    1341, 286, 243, 45};

constexpr double alpha = 0.2;    // 효율계수
constexpr double beta = 1.0 / 3; // 동원률

/*
	 * void simple_simulation(int data[M][K], doulbe (&time)[K])
	 *
	 * S개의 시나리오와 data로 부터 time을 산출
	 *
	 * (제약조건:
	 *   K가 무기의 종류지만, 우리는 0번째 무기만을
	 *   이용함. 오전 7시에 사건 발생, 하루의 일할 수
	 *   있는 시간은 8시간이고, 사건 발생 지점에 할당된
	 *   유회수기만을 사용하여 1시간 단위로 시뮬레이션
	 * )
	 * -----------------------------------------------
	 * Global Varaible
	 *
	 * int SCENARIO[S] : 유출유 사고 시나리오
	 *                   (index:지역, value: 유출유(톤))
	 * double alpha : 0.2 효율 계수
	 * double beta  : 1/3 동원률
	 *
	 * -----------------------------------------------
	 * Function Parameter
	 *
	 * int data[M][K]
	 *  : M개의 지점에 첫 번째 무기(유회수기)가 배치된
	 *    상태.
	 *
	 * (&time)[S]
	 *  : 시나리오로 부터 유회수기를 이용하여 유출유를
	 *    회수하는데 걸리는 시간을 산출
	 * -----------------------------------------------
	 *
	 * 예제 :
	 * // 첫 번째 지역에 사고 1톤,
	 * // 두 번째 지역에 사고 1톤,
	 * // ...
	 * // 네 번째 지역에 사고 4톤
	 * int SCENARIO[S] = { 1, 1, 2, 1 }
	 *
	 * // 첫 번째 지역에 유회수기 3개 배치,
	 * // ...
	 * // 네 번째 지역에 유회수기 1개 배치
	 * int data[4][1] = { {3},
	 *                    {6},
	 *                    {6},
	 *                    {1}
	 *                  }
	 *
	 * double time[4];
	 * simple_simulation(data, time);
	 *
	 * time == { 5, 3, 5, 31 }
	 *
	 */
void simple_simulation(const int data[M][K], double (&time)[S]);
void simulation(const int data[M][K], double (&time)[S], const int T_0 = 8, const int ship_speed = 10);

} // namespace sim
#endif /* SIMULATION_H */