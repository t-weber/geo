/**
 * @author Tobias Weber (orcid: 0000-0002-7230-1932)
 * @date may-2020
 * @license: see 'LICENSE' file
 */

#ifndef __CASTER_H__
#define __CASTER_H__

#include <QDialog>
#include <QTimer>
#include <QPaintEvent>
#include <QResizeEvent>
#include <QMouseEvent>
#include <QVector2D>

#include <memory>
#include <chrono>
#include <vector>
#include <array>

#include <boost/geometry.hpp>
namespace geo = boost::geometry;


#define NUM_RAYS 192


using t_real = float;
using t_vertex = geo::model::point<t_real, 2, geo::cs::cartesian>;
using t_poly = geo::model::polygon<t_vertex, true, false>;
using t_lines = geo::model::linestring<t_vertex>;


struct Casted
{
	t_real dist = -1.f;
	t_real column = -1.f;
	t_vertex vertex{};
};


class CasterWidget : public QWidget
{
public:
	using QWidget::QWidget;

	CasterWidget(QWidget *pParent);
	virtual ~CasterWidget();

	QPointF ToScreenCoords(const QVector2D& vec);
	QPointF ToSidescreenCoords(const QVector2D& vec);

protected:
	virtual void resizeEvent(QResizeEvent *pEvt) override;
	virtual void paintEvent(QPaintEvent *pEvt) override;
	virtual void mouseMoveEvent(QMouseEvent *pEvt) override;
	virtual void keyPressEvent(QKeyEvent* pEvt) override;
	virtual void keyReleaseEvent(QKeyEvent* pEvt) override;

	void tick(const std::chrono::milliseconds& ms);

private:
	t_real m_screenDims[2] = { -1, -1 };
	QTimer m_timer{};
	QPointF m_posMouse{};

	std::vector<t_poly> m_geo{};

	bool m_up{}, m_down{}, m_left{}, m_right{};

	// camera position and angle
	QVector2D m_pos{}, m_dir{};
	t_real m_angle{}, m_fov{};
	QVector2D m_fovlines[2];

	std::array<Casted, NUM_RAYS> m_casted{};

protected slots:
	void tick();
};


class CasterDlg : public QDialog
{
public:
	using QDialog::QDialog;
	CasterDlg(QWidget* pParent);
	~CasterDlg() = default;

private:
	std::shared_ptr<CasterWidget> m_pWidget;
};


#endif
