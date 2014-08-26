#include "Loader.h"

#include <dlib/xml_parser.h>
#include <iostream>
#include <fstream>

class doc_handler : public dlib::document_handler
{

public:
    virtual void start_document ();
    virtual void end_document ();

    virtual void start_element (const unsigned long line_number,
                                const std::string& name,
                                const dlib::attribute_list& atts);
    virtual void end_element (const unsigned long line_number,
                              const std::string& name);
    virtual void characters (const std::string& data);

    virtual void processing_instruction (const unsigned long line_number,
                                         const std::string& target,
                                         const std::string& data);
};

struct Loader::PImpl
{
public:
    doc_handler m_docHandler;
};

Loader::Loader()
{
    m_impl = new PImpl();
}

Loader::~Loader()
{
    delete m_impl;
}


void doc_handler::start_document()
{
    std::cout << "parsing begins" << std::endl;
}

void doc_handler::end_document ()
{
    std::cout << "Parsing done" << std::endl;
}

void doc_handler::start_element (const unsigned long line_number,
                                         const std::string& name,
                                         const dlib::attribute_list& atts)
{
    std::cout << "on line " << line_number << " we hit the <" << name << "> tag" << std::endl;

    // print all the tag's attributes
    atts.reset();
    while (atts.move_next())
    {
        std::cout << "\tattribute: " << atts.element().key() << " = " << atts.element().value() << std::endl;
    }
}

void doc_handler::end_element (const unsigned long line_number,
                                       const std::string& name)
{
    std::cout << "on line " << line_number << " we hit the closing tag </" << name << ">" << std::endl;
}

void doc_handler::characters (const std::string& data)
{
    std::cout << "Got some data between tags and it is:\n" << data << std::endl;
}

void doc_handler::processing_instruction (const unsigned long line_number,
                                     const std::string& target,
                                     const std::string& data)
{
    std::cout << "on line " << line_number << " we hit a processing instruction with a target of '"
        << target << "' and data '" << data << "'" << std::endl;
}
