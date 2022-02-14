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
#ifndef STATANALY_D_DOUBLE_DISPATCHER_H_
#define STATANALY_D_DOUBLE_DISPATCHER_H_

#include "type_info.h"
#include <map>


namespace statanaly {

/**
 * @brief Backend to FnDispatcher.
 * 
 * Reference: Modern C+ Design by Andrei Alexandrescu
 * 
 * @tparam BaseLhs 
 * @tparam BaseRhs 
 * @tparam ResultType 
 * @tparam (*)(BaseLhs&, BaseRhs&)
 * @see FnDispatcher 
 */
template <class BaseLhs,
		class BaseRhs = BaseLhs,
		typename ResultType = void,
		typename CallbackType = ResultType (*)(BaseLhs&, BaseRhs&)>
class BasicDispatcher {
	typedef std::pair<TypeInfo,TypeInfo> KeyType;
	typedef CallbackType MappedType;
	typedef std::map<KeyType, MappedType> MapType;

public:
	template <class SomeLhs, class SomeRhs>
	void add(CallbackType fun) {
		const KeyType key(typeid(SomeLhs), typeid(SomeRhs));
		callbackMap_[key] = fun;
	};

	/* Search and Invocation */
	ResultType go(BaseLhs& lhs, BaseRhs& rhs) {
		const KeyType key(typeid(lhs), typeid(rhs));
		auto i = callbackMap_.find(KeyType(typeid(lhs), typeid(rhs)));
		if (i == callbackMap_.end()) {
			throw std::runtime_error("Function not found");
		}

		return (i->second)(lhs, rhs);
	}

private:
	MapType callbackMap_;
};


/**
 * @brief Double Dispatcher.
 * 
 * Reference: Modern C+ Design by Andrei Alexandrescu
 * 
 * @tparam BaseLhs 
 * @tparam BaseRhs 
 * @tparam ResultType 
 */
template <class BaseLhs,
		class BaseRhs = BaseLhs,
		typename ResultType = void>
class FnDispatcher {
private:
	BasicDispatcher<BaseLhs, BaseRhs, ResultType> backEnd_;

public:
	template <class ConcreteLhs,
			class ConcreteRhs,
			ResultType (*callback)(ConcreteLhs&, ConcreteRhs&),
			bool symmetric=true>
	void add() {
		struct Local {
			/**
			 * @brief A trampoline function is saved in the lookup table as the callback.
			 * 
			 * The trampoline function has the concrete types information saved during registration.
			 */
			static ResultType Trampoline(BaseLhs& lhs, BaseRhs& rhs) {
				return callback(
						dynamic_cast<ConcreteLhs&>(lhs),
						dynamic_cast<ConcreteRhs&>(rhs));
			}
			// symmetry support
			static ResultType TrampolineR(BaseRhs& rhs, BaseLhs& lhs) {
				return Trampoline(lhs,rhs);
			}
		};

		backEnd_.template add<ConcreteLhs, ConcreteRhs>(&Local::Trampoline);
		
		if (symmetric) {
			// symmetry support
			backEnd_.template add<ConcreteRhs, ConcreteLhs>(&Local::TrampolineR);
		}
	}

	ResultType go(BaseLhs& lhs, BaseRhs& rhs) {
		return backEnd_.go(lhs,rhs);
	}
};

}   // namespace statanaly

#endif