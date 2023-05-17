#include <gtest/gtest.h>
#include <tinyxml2/tinyxml2.h>

namespace TestBlob
{
    namespace
    {
        TEST(TestBlob, tinyxml2_load)
        {
            double dblMass = 0.0;
            tinyxml2::XMLError xmlErr;
            tinyxml2::XMLDocument xmlDoc;
            xmlErr = xmlDoc.LoadFile("model/quadip/QuadIPConfig.xml");
            if (tinyxml2::XML_SUCCESS != xmlErr)
            {
                std::cout << "tinyxml2::Error (" << xmlErr << ")" << std::endl;
            }
            tinyxml2::XMLElement *xmlElemRoot = xmlDoc.RootElement();
            tinyxml2::XMLElement *elem = xmlElemRoot->FirstChildElement();
            for (; elem != NULL; elem = elem->NextSiblingElement())
            {
                if (!strcmp("mass", elem->Name()))
                {
                    dblMass = elem->DoubleAttribute("value");
                }
            }
            std::cout << "mass = " << dblMass << std::endl;
        }
    }
}