#ifndef __UNION_FIND__
#define __UNION_FIND__

#include <cstdint>
#include <vector>

struct union_find {
	std::vector< int64_t > e;
	union_find(int64_t n) { e.assign(n, -1); }
	int64_t find_set(int64_t x) {
		return (e[x] < 0 ? x : e[x] = find_set(e[x]));
	}
	bool same_set (int64_t x, int64_t y) { return find_set(x) == find_set(y); }
	int64_t size (int64_t x) { return -e[find_set(x)]; }
	bool union_set (int64_t x, int64_t y) {
		x = find_set(x), y = find_set(y);
		if (x == y) return 0;
		if (x > y) std::swap(x, y);
		e[x] += e[y], e[y] = x;
		return 1;
	}
	void clear() { e.clear(); }
};
#endif
