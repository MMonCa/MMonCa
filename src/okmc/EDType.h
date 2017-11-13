/*
 * EDType.h
 *
 *  Created on: Jun 29, 2011
 *
 * Author: ignacio.martin@imdea.org
 *
 * Copyright 2014 IMDEA Materials Institute, Getafe, Madrid, Spain
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *    http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 *      Virtual class to allow using user-specified inputs for extended defects or hard-coded definitions
 *      It contains two classes: EDType (virtual base class) and EDTypeGeneric (user specified values)
 */

#ifndef EDTYPE_H_
#define EDTYPE_H_

#include "kernel/ParticleType.h"
#include "kernel/Coordinates.h"
#include "kernel/RNG.h"
#include "kernel/IDContainer.h"
#include "io/ArrheniusAlloys.h"
#include <vector>
#include <set>

namespace OKMC
{
	class Particle;

	class EDType
	{
	public:
		void checkAxes();

		enum ASPECT_TYPE { aspect_disk, aspect_irregular, aspect_void, aspect_plane311 };
		enum MIGRATION   { mig_3d, mig_parallel, mig_perpendicular };

		virtual Kernel::Coordinates surface  (Kernel::RNG &, const Kernel::Coordinates c[], //creation point
				const Kernel::Coordinates &center, const Kernel::ID &, const Particle *) const = 0;
		virtual const Particle * emitFrom (Kernel::RNG &, const Kernel::Coordinates &center, //emission point.
				const std::vector<Particle *> &parts) const = 0;
		virtual ASPECT_TYPE      getAspect() const = 0;

		std::string _name;
		float _densityNm;
		unsigned _transformTo;
		unsigned _transformFrom;
		bool _ivmodel;

		Kernel::Coordinates _axes[3];
		Kernel::Coordinates _doNotCreateIn;
		MIGRATION _migration;
		float _lambda;

		//fuzz model
		Kernel::P_TYPE _expand_impurity;
		float          _expand_impurity_volume_nm3;
		float          _expand_capture_radius_P;  //distance to surface for fuzz = P + M*size
		float          _expand_capture_radius_M;

		bool _percolation;

		struct CLType
		{
			CLType();
			IO::ArrheniusAlloys _eForm;
			float _pref[Kernel::MAX_PARTICLES]; //one per possible particle emitted.
			bool _interactMP[Kernel::MAX_PARTICLES]; //interaction with particles. One by one
			bool _interactMP_sink[Kernel::MAX_PARTICLES]; //whether the MP is a sink.
			float _interactMP_radius[Kernel::MAX_PARTICLES];
			IO::ArrheniusAlloys _arr[4]; // 0->mig, 1->transform.to, 2->transform.from, 3->recomb
		};

		Kernel::IDContainer<CLType> _hash;
		std::map<unsigned, float> _interactMC; //interaction between clusters, the whole set.
			// the float is the capture distance between particles.
	};
}

#endif /* EDTYPE_H_ */
