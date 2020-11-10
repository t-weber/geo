/**
 * convex hull test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 15-Aug-2020
 * @license: see 'LICENSE' file
 */

#ifndef __VIS_GUI_H__
#define __VIS_GUI_H__


#include <QMainWindow>
#include <QLabel>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QGraphicsItem>

#include <memory>
#include <vector>


class Vertex : public QGraphicsItem
{
public:
	Vertex(const QPointF& pos, double rad = 15.);
	virtual ~Vertex();

	virtual QRectF boundingRect() const override;
	virtual void paint(QPainter*, const QStyleOptionGraphicsItem*, QWidget*) override;

private:
	double m_rad = 15.;
};


class VisView : public QGraphicsView
{Q_OBJECT
public:
	VisView(QGraphicsScene *scene=nullptr, QWidget *parent=nullptr);
	virtual ~VisView();

	VisView(VisView&) = delete;
	const VisView& operator=(const VisView&) const = delete;

	void AddVertex(const QPointF& pos);
	void ClearVertices();
	const std::vector<Vertex*>& GetVertices() const { return m_vertices; }

	void UpdateAll();
	void UpdateEdges();

protected:
	virtual void mousePressEvent(QMouseEvent *evt) override;
	virtual void mouseReleaseEvent(QMouseEvent *evt) override;
	virtual void mouseMoveEvent(QMouseEvent *evt) override;

	virtual void resizeEvent(QResizeEvent *evt) override;

private:
	QGraphicsScene *m_scene = nullptr;

	std::vector<Vertex*> m_vertices{};
	std::vector<QGraphicsItem*> m_edges{};

	bool m_dragging = false;

signals:
	void SignalMouseCoordinates(double x, double y);
};



class VisWnd : public QMainWindow
{
public:
	using QMainWindow::QMainWindow;

	VisWnd(QWidget* pParent = nullptr);
	~VisWnd();

	void SetStatusMessage(const QString& msg);

private:
	virtual void closeEvent(QCloseEvent *) override;

private:
	std::shared_ptr<QGraphicsScene> m_scene;
	std::shared_ptr<VisView> m_view;
	std::shared_ptr<QLabel> m_statusLabel;
};


#endif
