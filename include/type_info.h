/*
   Copyright 2022, Ansel Blumers

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
#ifndef STATANALY_D_TYPE_INFO_H_
#define STATANALY_D_TYPE_INFO_H_

#include <typeinfo>


namespace statanaly {

class TypeInfo {
public:
	TypeInfo() {
		// needed for containers
		pInfo_ = nullptr;
	};

	TypeInfo(const std::type_info& ti) {
		pInfo_ = &ti;
	};

	TypeInfo(const TypeInfo& fic) {
		pInfo_ = fic.pInfo_;
	};

	inline TypeInfo& operator=(const TypeInfo& fic) {
		this->pInfo_ = &fic.get();
		return *this;
	};

	inline bool before(const TypeInfo& rhs) const {
		return pInfo_->before(rhs.get());
	};

	inline const char* name() const {
		return pInfo_->name();
	};

	inline const std::type_info& get() const {
		return *pInfo_;
	}

private:
	const std::type_info* pInfo_;
};


// Need those comparators because TypeInfo is the key-type in std::map.
bool operator==(const TypeInfo& lhs, const TypeInfo& rhs);
bool operator!=(const TypeInfo& lhs, const TypeInfo& rhs);
bool operator< (const TypeInfo& lhs, const TypeInfo& rhs);
bool operator> (const TypeInfo& lhs, const TypeInfo& rhs);
bool operator<=(const TypeInfo& lhs, const TypeInfo& rhs);
bool operator>=(const TypeInfo& lhs, const TypeInfo& rhs);


}   // namespace statanaly

#endif