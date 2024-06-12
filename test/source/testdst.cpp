#include <doctest/doctest.h>
#include "dst/dst.hpp"
#include "dst/version.hpp"

#include <string>

TEST_CASE("DST") {
  using namespace dst;

  DST greeter("Tests");

  CHECK(greeter.greet(LanguageCode::EN) == "Hello, Tests!");
  CHECK(greeter.greet(LanguageCode::DE) == "Hallo Tests!");
  CHECK(greeter.greet(LanguageCode::ES) == "¡Hola Tests!");
  CHECK(greeter.greet(LanguageCode::FR) == "Bonjour Tests!");
}

TEST_CASE("DST version") {
  CHECK(std::string_view(DST_VERSION) == std::string_view("1.0"));
  CHECK(std::string(DST_VERSION) == std::string("1.0"));
}
