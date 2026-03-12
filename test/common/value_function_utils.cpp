#include <gtest/gtest.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include "util/value_function_utils.hpp"

namespace hcpwa_util = hcpwa::util;

TEST(value_function_utils, setValueFunction_and_getValueFunction) {
    hcpwa_util::ValueFunction vf;
    std::vector<double> data = {1.0, 2.0, 3.0};
    vf.set(0, 1, 2, 3, data, 3);
    std::vector<double> got = vf.get(0, 1, 2, 3);
    ASSERT_EQ(got.size(), 3u);
    EXPECT_DOUBLE_EQ(got[0], 1.0);
    EXPECT_DOUBLE_EQ(got[1], 2.0);
    EXPECT_DOUBLE_EQ(got[2], 3.0);
}

TEST(value_function_utils, setValueFunction_throws_on_wrong_size) {
    hcpwa_util::ValueFunction vf;
    std::vector<double> data = {1.0, 2.0};
    EXPECT_THROW(vf.set(0, 0, 0, 0, data, 3), std::invalid_argument);
}

TEST(value_function_utils, setValueFunction_throws_when_location_already_used) {
    hcpwa_util::ValueFunction vf;
    std::vector<double> data = {1.0, 2.0};
    vf.set(0, 0, 0, 0, data, 2);
    EXPECT_THROW(
        vf.set(0, 0, 0, 0, std::vector<double>{3.0, 4.0}, 2),
        std::runtime_error);
}

TEST(value_function_utils, getValueFunction_throws_when_not_found) {
    hcpwa_util::ValueFunction vf;
    EXPECT_THROW(vf.get(99, 99, 99, 99), std::invalid_argument);
}

TEST(value_function_utils, dumpValueFunctionToJson_writes_valid_structure) {
    hcpwa_util::ValueFunction vf;
    vf.set(0, 1, 2, 3, {1.0, 2.0}, 2);
    vf.set(0, 1, 2, 4, {3.0, 4.0, 5.0}, 3);
    vf.set(1, 0, 0, 0, {0.5}, 1);

    std::filesystem::path tmp = std::filesystem::temp_directory_path() / "value_function_utils_test.json";
    ASSERT_NO_THROW(vf.dumpToJson(tmp.string()));

    ASSERT_TRUE(std::filesystem::exists(tmp)) << "Output file was not created";

    std::ifstream f(tmp.string());
    std::ostringstream contents;
    contents << f.rdbuf();
    std::string json = contents.str();

    EXPECT_GT(json.size(), 0u);
    EXPECT_EQ(json.front(), '{');
    EXPECT_EQ(json.back(), '}');
    EXPECT_NE(json.find("\"0\":"), std::string::npos);
    EXPECT_NE(json.find("\"1\":"), std::string::npos);
    EXPECT_NE(json.find("[1,2]"), std::string::npos);
    EXPECT_NE(json.find("[3,4,5]"), std::string::npos);
    EXPECT_NE(json.find("[0.5]"), std::string::npos);

    std::filesystem::remove(tmp);
}

TEST(value_function_utils, dumpValueFunctionToJson_empty_map) {
    hcpwa_util::ValueFunction vf;
    std::filesystem::path tmp = std::filesystem::temp_directory_path() / "value_function_utils_empty_test.json";
    ASSERT_NO_THROW(vf.dumpToJson(tmp.string()));

    std::ifstream f(tmp.string());
    std::ostringstream contents;
    contents << f.rdbuf();
    std::string json = contents.str();
    EXPECT_EQ(json, "{}");
    std::filesystem::remove(tmp);
}

TEST(value_function_utils, dumpValueFunctionToJson_throws_on_bad_path) {
    hcpwa_util::ValueFunction vf;
    vf.set(0, 0, 0, 0, {1.0}, 1);
    EXPECT_THROW(vf.dumpToJson("/nonexistent/dir/out.json"),
                 std::runtime_error);
}

TEST(value_function_utils, dumpVectorToJson_writes_array) {
    std::vector<double> vec = {1.0, 2.5, 3.0};
    std::filesystem::path tmp =
        std::filesystem::temp_directory_path() / "dump_vector_test.json";
    ASSERT_NO_THROW(hcpwa_util::dumpVectorToJson(vec, tmp.string()));
    ASSERT_TRUE(std::filesystem::exists(tmp));

    std::ifstream f(tmp.string());
    std::ostringstream contents;
    contents << f.rdbuf();
    std::string json = contents.str();
    EXPECT_EQ(json.front(), '[');
    EXPECT_EQ(json.back(), ']');
    EXPECT_NE(json.find("2.5"), std::string::npos);
    EXPECT_GE(json.size(), 7u);  // at least "[1,2.5,3]" or similar
    std::filesystem::remove(tmp);
}

TEST(value_function_utils, dumpVectorToJson_empty_vector) {
    std::vector<double> vec;
    std::filesystem::path tmp =
        std::filesystem::temp_directory_path() / "dump_vector_empty_test.json";
    ASSERT_NO_THROW(hcpwa_util::dumpVectorToJson(vec, tmp.string()));
    std::ifstream f(tmp.string());
    std::ostringstream contents;
    contents << f.rdbuf();
    EXPECT_EQ(contents.str(), "[]");
    std::filesystem::remove(tmp);
}

TEST(value_function_utils, dumpVectorToJson_throws_on_bad_path) {
    std::vector<double> vec = {1.0};
    EXPECT_THROW(
        hcpwa_util::dumpVectorToJson(vec, "/nonexistent/dir/vector.json"),
        std::runtime_error);
}
