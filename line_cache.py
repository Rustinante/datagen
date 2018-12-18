class LineCache:
    def __init__(self):
        self.capacity = 1000
        self.count_per_eviction = 200
        self.cache = {}
        self.key_list = []

    def evict(self):
        for key in self.key_list[:self.count_per_eviction]:
            del self.cache[key]
        self.key_list = self.key_list[self.count_per_eviction:]

    def __contains__(self, item):
        return item in self.cache

    def __getitem__(self, item):
        return self.cache[item]

    def __setitem__(self, key, value):
        if len(self.cache) >= self.capacity:
            self.evict()
        self.cache[key] = value
        if key not in self.key_list:
            self.key_list.append(key)

    def __len__(self):
        return len(self.cache)

    def keys(self):
        return self.cache.keys()

    def values(self):
        return self.cache.values()

    def items(self):
        return self.cache.items()
