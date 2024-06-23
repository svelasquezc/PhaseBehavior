#ifndef GERG_TYPES_HPP
#define GERG_TYPES_HPP

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <vector>
#include <functional>
#include <string_view>

namespace PhaseBehavior::EoS::GERG::Types{
    constexpr std::pair<std::string_view, std::string_view> pairSort(std::pair<std::string_view, std::string_view> pair) {
        auto [first, second] = pair;
        if (first > second) std::swap(first, second);

        return {first, second};
    }
    
    // Custom hash function for a pair of strings
    struct PairHash {
        std::size_t operator()(std::pair<std::string_view, std::string_view> keyPair) const {
            // Create a vector of the two strings and sort it
            auto [firstKey, secondKey] = pairSort(keyPair);

            // Combine the hashes of the sorted keys
            return std::hash<std::string_view>()(firstKey) ^ std::hash<std::string_view>()(secondKey);
        }
    };

    // Custom equality function for a pair of strings
    struct PairEqual {
        bool operator()(const std::pair<std::string_view, std::string_view>& a, const std::pair<std::string_view, std::string_view>& b) const {
            // Create sorted vectors of the pairs
            auto sorted_a = pairSort(a);
            auto sorted_b = pairSort(b);

            // Compare the sorted vectors
            return sorted_a == sorted_b;
        }
    };

    template <typename ContainedType>
    using BinaryParameter = std::unordered_map<std::pair<std::string_view, std::string_view>, ContainedType, PairHash, PairEqual>;
}

#endif /* GERG_TYPES_HPP */