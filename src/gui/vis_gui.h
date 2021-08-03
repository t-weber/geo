/**
 * vis test program
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date 11-Nov-2020
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

#include "geo_algos.h"
using t_real = double;


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
	using t_vec = m::vec<t_real, std::vector>;
	using t_mat = m::mat<t_real, std::vector>;

	static const constexpr t_real g_eps = 1e-5;


public:
	VisView(QGraphicsScene *scene=nullptr, QWidget *parent=nullptr);
	virtual ~VisView();

	VisView(VisView&) = delete;
	const VisView& operator=(const VisView&) const = delete;

	void AddVertex(const QPointF& pos);
	void ClearVertices();
	const std::vector<Vertex*>& GetVertexElems() const { return m_elems_vertices; }

	void UpdateAll();
	void UpdateEdges();
	void UpdateSplitPolygon();
	void UpdateKer();

	void SetSortVertices(bool b);
	bool GetSortVertices() const { return m_sortvertices; }

	void SetCalcSplitPolygon(bool b);
	bool GetCalcSplitPolygon() const { return m_splitpolygon; }

	void SetCalcKernel(bool b);
	bool GetCalcKernel() const { return m_calckernel; }

protected:
	virtual void mousePressEvent(QMouseEvent *evt) override;
	virtual void mouseReleaseEvent(QMouseEvent *evt) override;
	virtual void mouseMoveEvent(QMouseEvent *evt) override;

	virtual void resizeEvent(QResizeEvent *evt) override;

private:
	QGraphicsScene *m_scene = nullptr;

	std::vector<Vertex*> m_elems_vertices{};
	std::vector<QGraphicsItem*> m_elems_edges{}, m_elems_ker{}, m_elems_split{};

	bool m_dragging = false;

	std::vector<t_vec> m_vertices{};

	bool m_sortvertices = true;
	bool m_splitpolygon = true;
	bool m_calckernel = true;

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
