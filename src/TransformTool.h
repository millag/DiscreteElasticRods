#ifndef TRANSFORMTOOL_H
#define TRANSFORMTOOL_H

#include "Camera.h"


class TransformHandle
{
public:

    enum TransformMode
    {
        TM_None = 0,
        TM_Rotation,
        TM_Scale,
        TM_Translation,
        TM_Count
    };

    TransformHandle(Camera* cam = nullptr):
        m_cam(cam)
    { }
    ~TransformHandle()
    { }

    bool isActive() const { return m_cam && m_mode != TM_None; }
    TransformMode getMode() const { return m_mode; }
    void setMode(TransformMode mode) { m_mode = mode; }

    void update(mg::Real dx, mg::Real dy)
    {
        if (!m_cam)
        {
            return;
        }
        switch (m_mode)
        {
            case TM_Rotation:
            {
                m_cam->rotate(dx, dy);
                break;
            }
            case TM_Translation:
            {
                m_cam->pan(dx, dy);
                break;
            }
            case TM_Scale:
            {
                m_cam->zoom(dx, dy);
                break;
            }
        }
    }

private:
/// @brief transformation mode
    TransformMode m_mode = TM_None;
/// @brief transform object
    Camera* m_cam = nullptr;
};

class QObject;
class INode;
class IViewport;
class IViewportRenderer;

enum class UIEventType
{
    UIEventUndefined,
    UIEventMousePress,
    UIEventMouseDrag,
    UIEventMouseRelease,
    UIEventMouseClick,
    UIEventMouseDoubleClick,
    UIEventHoverEnter,
    UIEventHoverLeave,
    UIEventKeyPress,
    UIEventKeyRelease,
    UIEventKeyPressed,
    UIEventTouch,
    UIEventWheel,
    UIEventGesture,

    UIEvent_Count,
};


class IInteraction
{
public:
    virtual ~IInteraction()
    { }

    /// interaction classification type such as SELECT, TRANSFORM, CREATE, MODIFY, etc.
    virtual int getType() const = 0;
    /// unique id of the interation plugin
    virtual int getID() const = 0;

    /// descriptive name, for debugging purposes
    virtual const char* getDescriptiveName() const { return nullptr; }

    /// called when the interation should become active/inactive
    virtual bool isActive() const = 0;
    virtual int activate( IInteraction* lastInteraction ) = 0;
    virtual void deactivate( IInteraction* nextInteraction ) = 0;

    /// make DG Node associations
    virtual void attachNode( INode* node ) = 0;
    virtual void detachNode( INode* node ) = 0;

    /// register user interaction callbacks
    /// use IViewport::registerUICallback() methods to register callbacks
    /// for user input events that the Interaction is interested in
    /// callback must satisfy UIInterationCallback signature
    virtual void registerUICallbacks( IViewport& ) = 0;
    virtual void deregisterUICallbacks( IViewport& ) = 0;

    /// viewport drawing
    virtual bool supportsGizmo( IViewport& viewport ) const = 0;
    virtual void preRenderGizmo( IViewport& viewport ) = 0;
    virtual void renderGizmo( IViewportRenderer& renderer ) = 0;
};


class BarBase
{
    QObject* m_source;
    QObject* m_receiver;
};

class Bar
{
public:
    UIEventType getType() const { return m_type; }
    const QObject* getQtEvent() const { return m_qtevent; }
    uint64_t getTimestamp() const { return m_timestamp; }

    void example(BarBase& base)
    {
        getType();
    }

    inline static bool isUIEvent( UIEventType type )
    {
        return type == UIEventType::UIEventMousePress;
    }

private:
    Bar() = default;
    ~Bar() = default;

    UIEventType m_type;
    QObject* m_qtevent;
    uint64_t m_timestamp;

    friend class FooImplPrivate;
};

class Foo
{
public:
    static Foo* create();
    static void destroy(Foo* foo);

    bool isValid() const { return m_pimpl != nullptr; }
    void destroy();

private:
    Foo() = default;
    ~Foo() = default;

    struct FooPImpl;
    FooPImpl* m_pimpl;
    friend class FooImplPrivate;
};

#endif // TRANSFORMTOOL_H
