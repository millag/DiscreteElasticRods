#include "TransformTool.h"
#include <iostream>

struct Foo::FooPImpl
{
    FooPImpl() { }
    ~FooPImpl() { }

    std::vector<mg::Matrix4D> m_tm;
};


Foo* Foo::create()
{
    auto foo = new Foo();
    foo->m_pimpl = new FooPImpl();
    return foo;
}

void Foo::destroy(Foo* foo)
{
    if (foo)
    {
        foo->destroy();
        delete foo;
    }
}

void Foo::destroy()
{
    if (m_pimpl)
    {
        delete m_pimpl;
        m_pimpl = nullptr;
    }
}
