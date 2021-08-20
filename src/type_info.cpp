#include "type_info.h"


namespace statanaly {

// Need those comparators because TypeInfo is the key-type in std::map.
bool operator==(const TypeInfo& lhs, const TypeInfo& rhs) {
	return lhs.get() == rhs.get();
};
bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs) {
	return !(lhs == rhs);
};
bool operator< (const TypeInfo& lhs, const TypeInfo& rhs) {
	return lhs.before(rhs);
};
bool operator> (const TypeInfo& lhs, const TypeInfo& rhs) {
	return rhs < lhs;
};
bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs) {
	return !(lhs > rhs);
};
bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs) {
	return !(lhs < rhs);
};


} 	// namespace statanaly