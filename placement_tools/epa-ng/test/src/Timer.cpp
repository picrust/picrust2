#include "Epatest.hpp"

#include <chrono>
#include <thread>

#include "util/Timer.hpp"

TEST(Timer, pause)
{
  Timer<std::chrono::milliseconds> t;
  t.start();

  std::this_thread::sleep_for(std::chrono::milliseconds(2));

  t.pause();

  std::this_thread::sleep_for(std::chrono::milliseconds(1));

  t.resume();
  t.stop();

  ASSERT_DOUBLE_EQ(t.average(), 2.0);

}

TEST(Timer, construct_from_avg)
{
  Timer<> t;
  t.start();

  std::this_thread::sleep_for(std::chrono::milliseconds(2));

  t.stop();

  Timer<> tt(t.avg_duration());

  ASSERT_DOUBLE_EQ(tt.average(), t.average());

}