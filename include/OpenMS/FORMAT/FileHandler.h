// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_FILEHANDLER_H
#define OPENMS_FORMAT_FILEHANDLER_H

#include <OpenMS/config.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/FORMAT/DTA2DFile.h>
#include <OpenMS/FORMAT/MzXMLFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/FORMAT/ANDIFile.h>

namespace OpenMS
{
	/**
		 @brief Facilitates file handling by file type recognition.

		 This class provides file type recognition from the file name and
		 from the file content.

		 It also offer a common interface to load MSExperiment data
		 and allows querying for supported file types.

		 @ingroup FileIO
	*/
	class FileHandler
	{
	 public:
		/**
			 @brief The known file types.

			 @see nameToType and typeToName
		*/
		enum Type
			{
				UNKNOWN,        ///< Unknown file extension
				DTA,            ///< DTA file (.dta)
				DTA2D,          ///< DTA2D file (.dta2d)
				MZDATA,         ///< MzData file (.MzData)
				MZXML,          ///< MzXML file (.MzXML)
				FEATURE,        ///< OpenMS feature file (.feat)
				ANDIMS,         ///< ANDI\\MS file (.cdf)
				SIZE_OF_TYPE    ///< No file type. Simply stores the number of types
			};

		/// String representations of the file types
		static const std::string NamesOfTypes[SIZE_OF_TYPE];

		/// Determines the file type from a file name
		Type getTypeByFileName(const String& filename);

		/// Determines the file type of XML files by parsing the first few lines
		Type getTypeByContent(const String& filename) throw (Exception::FileNotFound);

		/// Converts a file type name into a Type
		Type nameToType(const String& name);

		/// Converts a Type into a file type name
		String typeToName(Type type);

		/// Returns if the file type is supported in this build of the library
		bool isSupported(Type type);

		/**
			 @brief Loads a file into an MSExperiment

			 @param filename the Filename of the file to load.
			 @param exp The MSExperiment to load the data into.
			 @param force_type Forces to load the file with that file type.<BR>
			 If no type is forced, it is determined from the extention ( or from the content if that fails).

			 @return true if the file could be loaded, false otherwise
		*/
		template <class PeakType> bool loadExperiment(const String& filename, MSExperiment<PeakType>& exp, Type force_type = UNKNOWN)
		{
			Type type;
			if (force_type != UNKNOWN)
			{
				type = force_type;
			}
			else
			{
				type = getTypeByFileName(filename);
				try
				{
					if (type == UNKNOWN)
					{
						type = getTypeByContent(filename);
					}
				}
				catch(Exception::FileNotFound)
				{
					// ???? I assume this is the intended behavior?  (Clemens asking Marc, 2006-01-11)
					return false;
				}
			}

			//load right file
			switch(type)
			{
			case DTA:
				exp.reset();
				exp.resize(1);
				DTAFile().load(filename,exp[0]);
				return true;
				break;
			case DTA2D:
				DTA2DFile().load(filename,exp);
				return true;
				break;
			case MZXML:
				MzXMLFile().load(filename,exp);
				return true;
				break;
			case MZDATA:
				MzDataFile().load(filename,exp);
				return true;
				break;
#ifdef ANDIMS_DEF
			case ANDIMS:
				ANDIFile().load(filename,exp);
				return true;
				break;
#endif
			default:
				return false;
			}
		}
	};

} //namespace

#endif //OPENMS_FORMAT_FILEHANDLER_H
