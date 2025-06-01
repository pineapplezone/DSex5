// by 11227205 資訊二乙 劉至嘉 & 11027214 楊碕萍.
#include <algorithm>
#include <cctype>
#include <charconv>
#include <concepts>
#include <expected>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <ios>
#include <iostream>
#include <map>
#include <optional>
#include <print>
#include <queue>
#include <ranges>
#include <stack>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
    // inspired (copied) from absl::StatusCode
    // go to https://abseil.io/docs/cpp/guides/status-codes for documentation
    // NOLINTBEGIN
    enum struct StatusCode : int {
        kOk = 0,
        kCancelled = 1,
        kUnknown = 2,
        kInvalidArgument = 3,
        kDeadlineExceeded = 4,
        kNotFound = 5,
        kAlreadyExists = 6,
        kPermissionDenied = 7,
        kResourceExhausted = 8,
        kFailedPrecondition = 9,
        kAborted = 10,
        kOutOfRange = 11,
        kUnimplemented = 12,
        kInternal = 13,
        kUnavailable = 14,
        kDataLoss = 15,
        kUnauthenticated = 16,
    };
    // NOLINTEND

    template <typename T>
        requires std::integral<T> || std::floating_point<T>
    constexpr T StrTo(std::string_view str) noexcept {
        T value{};
        std::from_chars(str.data(), str.data() + str.size(), value);
        return value;
    }

    template <typename T>
        requires std::integral<T> || std::floating_point<T>
    constexpr std::optional<T> TryStrTo(std::string_view str) noexcept {
        T value{};
        const auto [ptr, errc] =
            std::from_chars(str.data(), str.data() + str.size(), value);
        if (errc == std::errc{}) {
            return value;
        }
        return {};
    }

    template <typename T>
    std::optional<T> Scan(std::istream& in = std::cin) noexcept {
        T val;
        in >> val;
        if (in.fail()) [[unlikely]] {
            return {};
        }
        return val;
    }

    template <typename T>
    std::optional<T> Scan(std::string_view prompt) {
        std::print("{}", prompt);
        return Scan<T>();
    }

    std::string ScanString(std::string_view prompt) {
        std::print("{}", prompt);
        std::string val;
        std::cin >> val;
        return val;
    }

    struct Message {
        std::string id1;
        std::string id2;
        float weight{};
    };

    class AdjacencyList {
    public:
        struct Edge {
            auto operator<=>(const Edge& other) const { return id <=> other.id; }
            auto operator==(const Edge& other) const { return id == other.id; }

            std::string_view id;
            float weight{};
        };

        using VertexType = std::string_view;
        using EdgesType = std::vector<Edge>;

        using ValueType = std::map<VertexType, EdgesType>;
        using MappedType = std::pair<VertexType, EdgesType>;

        using ConstIterator = ValueType::const_iterator;
        using Iterator = ValueType::iterator;

        void Insert(const Message& val) {
            const auto& [id1, id2, weight] = val;
            if (data_.contains(id1)) {
                data_[id1].emplace_back(id2, weight);
                std::ranges::sort(data_[id1]);
            }
            else {
                data_.emplace(id1, EdgesType{ {.id = id2, .weight = weight} });
            }

            if (data_.contains(id2)) {
                data_[id2].emplace_back(id1, weight);
                std::ranges::sort(data_[id2]);
            }
            else {
                data_.emplace(id2, EdgesType{ {.id = id1, .weight = weight} });
            }
        }
        void Clear() { data_.clear(); }

        float edge_weight(const std::string_view from,
            const std::string_view to) const {
            for (const auto& vec = data_.at(from); const auto [id, weight] : vec) {
                if (id == to) {
                    return weight;
                }
            }
        }

        bool Contains(const std::string_view id) const { return data_.contains(id); }
        bool IsEmpty() const { return data_.empty(); }
        size_t size() const { return data_.size(); }
        size_t nodes() const {
            size_t nodes{};
            for (const auto& val : data_ | std::views::values) {
                nodes += val.size();
            }
            return nodes;
        }
        Iterator begin() { return data_.begin(); }
        Iterator end() { return data_.end(); }
        ConstIterator begin() const { return data_.cbegin(); }
        ConstIterator end() const { return data_.cend(); }

        const EdgesType& at(const VertexType& id) const { return data_.at(id); }
        EdgesType& at(const VertexType& id) { return data_.at(id); }

    private:
        ValueType data_;
    };

    class ConnectedComponentAnalyzer {
    public:
        using Component = std::vector<std::string_view>;

        void Compute(const AdjacencyList& graph) {
            components_.clear();
            std::unordered_set<std::string_view> visited;

            for (const auto& id : graph | std::views::keys) {
                if (!visited.contains(id)) {
                    Component component;
                    DFS(id, graph, visited, component);
                    std::ranges::sort(component);
                    components_.push_back(std::move(component));
                }
            }

            std::ranges::sort(components_, [](const Component& a, const Component& b) {
                if (a.size() != b.size()) return a.size() > b.size();
                return a.front() > b.front();
                });
        }

        void PrintResults(std::ostream& out) const {
            std::println("<<< There are {} connected components in total. >>>",
                components_.size());

            std::println(out, "<<< There are {} connected components in total. >>>",
                components_.size());

            for (std::size_t count = 1; const auto & component : components_) {
                std::println("{{{:>2}}} Connected Component: size = {}", count,
                    component.size());

                std::println(out, "{{{:>2}}} Connected Component: size = {}", count,
                    component.size());
                for (std::size_t i = 0; i < component.size(); ++i) {
                    std::print(out, " \t({:>3}) {}", i + 1, component[i]);
                    if ((i + 1) % 8 == 0) {
                        std::println(out);
                    }
                }
                std::println(out);
                ++count;
            }

            std::cout << "\n";
        }

        Component& Find(const std::string_view val) {
            for (auto& component : components_) {
                if (std::ranges::contains(component, val)) {
                    return component;
                }
            }
        }

    private:
        static void DFS(std::string_view id, const AdjacencyList& graph,
            std::unordered_set<std::string_view>& visited,
            Component& component) {
            visited.insert(id);
            component.emplace_back(id);

            for (const auto [neighbor, weight] : graph.at(id)) {
                if (!visited.contains(neighbor)) {
                    DFS(neighbor, graph, visited, component);
                }
            }
        }

        std::vector<Component> components_;
    };

    class GraphSystem {
    public:
        static constexpr std::string_view kPrompt =
            "**********  Graph data applications  *********\n"
            "* 1. Build a graph and connected components  *\n"
            "* 2. Find shortest paths by Dijkstra         *\n"
            "**********************************************\n"
            "Input a choice(0, 1, 2) [0: QUIT]: ";

        StatusCode ExecuteCommand(const std::string_view input) {
            if (!IsValidCommand(input)) {
                std::println("\nCommand does not exist!\n");
                return StatusCode::kInvalidArgument;
            }
            using enum StatusCode;
            switch (const auto command = StrTo<int>(input); command) {
            case 0:
                return kAborted;
            case 1:
                BuildGraph();
                CalculateConnectedComponents();
                return kOk;
            case 2:
                FindShortestPath();
                return kOk;
            default:
                std::println("\nCommand does not exist!\n");
                return kOutOfRange;
            }
        }

    private:
        static bool IsValidCommand(std::string_view input) {
            return std::ranges::all_of(input, isdigit);
        }

        static std::vector<Message> ReadRelations(const std::string& file_name,
            const float threshold) {
            std::ifstream in_file{ file_name, std::ios::binary };

            char id_str[12]{};
            Message temp;
            std::vector<Message> relations;

            while (in_file.read(reinterpret_cast<char*>(&id_str), sizeof(id_str))) {
                temp.id1 = id_str;
                in_file.read(reinterpret_cast<char*>(&id_str), sizeof(id_str));
                temp.id2 = id_str;
                in_file.read(reinterpret_cast<char*>(&temp.weight), sizeof(temp.weight));

                if (temp.weight <= threshold) {
                    relations.push_back(std::move(temp));
                }
            }
            return relations;
        }

        static float ScanRealNumber() {
            auto input_number = ScanString("\nInput a real number in (0,1]: ");
            auto real_number = TryStrTo<float>(input_number);

            auto out_of_range = [](const float i) -> bool { return i <= 0 || i > 1; };
            while (!real_number || out_of_range(*real_number)) {
                if (real_number && out_of_range(*real_number) && *real_number >= 0) {
                    std::println("\n### It is NOT in (0,1] ###");
                }
                input_number = ScanString("\nInput a real number in (0,1]: ");
                real_number = TryStrTo<float>(input_number);
            }
            return *real_number;
        }

        void BuildGraph() {
            messages_.clear();
            graph_.Clear();

            const float real_number = ScanRealNumber();
            const auto file_number = ScanString("\nInput a file number ([0] Quit): ");

            if (const auto file_name = std::format("pairs{}.bin", file_number);
                std::filesystem::exists(file_name) && file_number != "0") {
                messages_ = ReadRelations(file_name, real_number);
                for (const auto& message : messages_) {
                    graph_.Insert(message);
                }
                const auto nodes = graph_.nodes();
                std::println(
                    "\n<<< There are {} IDs in total. >>>\n\n<<< There are {} nodes in "
                    "adjacency lists. >>>\n",
                    graph_.size(), nodes);

                if (real_number == 1) {
                    std::ofstream out{ std::format("pairs{}_1..adj", file_number) };
                    PrintList(out, nodes);
                }
                else {
                    std::ofstream out{
                        std::format("pairs{}_{}.adj", file_number, real_number) };
                    PrintList(out, nodes);
                }
            }

            file_number_ = file_number;
            threshold_ = real_number;
        }

        void CalculateConnectedComponents() {
            if (graph_.IsEmpty()) {
                std::println("### There is no graph and try it again. ###\n");
                return;
            }
            analyzer_.Compute(graph_);

            if (threshold_ == 1) {
                std::ofstream out{ std::format("pairs{}_1..cc", file_number_) };
                analyzer_.PrintResults(out);
            }
            else {
                std::ofstream out{ std::format("pairs{}_{}.cc", file_number_, threshold_) };
                analyzer_.PrintResults(out);
            }
        }
        
        static void WritePath(std::ostream& out, std::string_view origin,
            std::span<const AdjacencyList::Edge> data) {
            std::println(out, "\norigin: {}", origin);
            for (std::size_t count = 1; const auto [id, weight] : data) {
                std::print(out, "({:>2}) \t{}, {}\t", count, id, weight);

                if (count % 8 == 0) {
                    std::println(out);
                }
                ++count;
            }
            std::println(out);
        }

        void FindShortestPath() {
            if (graph_.IsEmpty()) {
                std::println("### There is no graph and try it again. ###\n");
                return;
            }
            std::ofstream out{ std::format("pairs{}_{}.ds", file_number_, threshold_),
                              std::ios::ate };
            auto id = ReadCorrectId();

            while (true) {
                if (id.has_value()) {
                    const auto path = DijkstraShortestPath(*id);
                    WritePath(out, *id, path);
                }
                else {
                    switch (id.error()) {
                        using enum StatusCode;

                    case kCancelled:
                        return;
                    case kNotFound:
                        std::println("\n### the student id does not exist! ###\n");
                    default:
                        break;
                    }
                }
                id = ReadCorrectId();
            }

            std::println();
        }

        static constexpr float kInfinite = std::numeric_limits<float>::max();

        std::vector<AdjacencyList::Edge> DijkstraShortestPath(std::string_view start) {
            using Edge = AdjacencyList::Edge;
            using Entry = std::pair<float, std::string_view>; // (cost, id)

            std::unordered_map<std::string_view, float> dist;
            std::unordered_set<std::string_view> visited;
            std::priority_queue<Entry, std::vector<Entry>, std::greater<>> pq;
            std::vector<Edge> result;

            // 初始：所有點設為無限大，起點設為0
            for (auto id : analyzer_.Find(start)) {
                dist[id] = kInfinite;
            }
            dist[start] = 0.0f;
            pq.emplace(0.0f, start);

            while (!pq.empty()) {
                auto [cur_cost, cur] = pq.top(); pq.pop();
                if (visited.contains(cur)) continue;
                visited.insert(cur);

                // 儲存目前節點與其距離
                if (cur != start) result.emplace_back(Edge{ .id = cur, .weight = cur_cost });

                for (const auto& [neighbor, weight] : graph_.at(cur)) {
                    if (!dist.contains(neighbor)) continue; // 不是同一個 component 的
                    if (visited.contains(neighbor)) continue;

                    float new_cost = cur_cost + weight;
                    if (new_cost < dist[neighbor]) {
                        dist[neighbor] = new_cost;
                        pq.emplace(new_cost, neighbor);
                    }
                }
            }

            // 按距離排序
            std::ranges::sort(result, [](const Edge& a, const Edge& b) {
                return a.weight < b.weight;
                });

            return result;
        }


        void PrintList(std::ostream& out, size_t nodes) const {
            std::println(out, "<<< There are {} IDs in total. >>>", graph_.size());
            for (std::size_t i = 1; const auto & v : graph_) {
                const auto& [id1, vec] = v;
                std::println(out, "[{:>3}] {}: ", i, id1);

                for (std::size_t j = 1; const auto & [id2, weight] : vec) {
                    std::print(out, "\t({:>2}) {},{:>7}", j, id2, weight);
                    if (j % 12 == 0) {
                        std::println(out);
                    }
                    ++j;
                }
                std::println(out);
                ++i;
            }
            std::println(out, "<<< There are {} nodes in adjacency lists. >>>", nodes);
        }

        std::expected<std::string, StatusCode> ReadCorrectId() const {
            using enum StatusCode;
            PrintAll();
            auto id = ScanString("\nInput a student ID [0: exit] ");
            if (id == "0") {
                return std::unexpected{ kCancelled };
            }
            if (!graph_.Contains(id)) {
                return std::unexpected{ kNotFound };
            }

            return id;
        }

        void PrintAll() const {
            for (const auto& id : graph_ | std::views::keys) {
                std::print("{:>12}", id);
            }
        }

        ConnectedComponentAnalyzer analyzer_;
        std::vector<Message> messages_;
        AdjacencyList graph_;
        std::string file_number_;
        float threshold_ = 0.0;
    };

}  // namespace

int main() {
    GraphSystem system;
    for (auto status = StatusCode::kOk; status != StatusCode::kAborted;) {
        const auto command_string = ScanString(GraphSystem::kPrompt);
        status = system.ExecuteCommand(command_string);
    }
}